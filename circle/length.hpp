// Observed-only, fixed-cutoff acceleration. The old interval kernel is unchanged.
// phi(q,r) = F_n(r) + eta ||ell(q)-ell(r)|| is jointly convex and
// 2-strongly convex in r. Every decision below is an outward interval proof.
namespace convex_length {
using namespace circle_kernel;
V root(V x) {
    if (x.hi < 0) return fc::whole();
    return {x.lo <= 0 ? 0 : fc::down(std::sqrt(x.lo)),
            fc::up(std::sqrt(std::max(0.0,x.hi)))};
}
double middle(V x) { return x.lo*.5+x.hi*.5; }
V risk(const Arc& a,V q) { return ad(sq(su(q,a.mean)),a.variance); }
struct Evaluation { V value, gr, gq, variance; bool differentiable=false; };
struct Engine {
    Data& data;
    std::vector<double> risk_lower, risk_witness;
    std::vector<std::vector<double>> warm;
    double fstar_upper=INFINITY;
    uint64_t evaluations=0, pairs=0, risk_skips=0, risk_out=0,
             diagonal_out=0, accepted=0, rejected=0, unresolved=0;
    Engine(Data& d):data(d) {
        for (const Arc& a:data.arcs) {
            V outer{a.left.lo,a.right.hi};
            double lower=ad(a.variance,sq(su(a.mean,outer))).lo;
            risk_lower.push_back(std::max(0.0,lower));
            double r=std::max(a.left.hi,std::min(a.right.lo,middle(a.mean)));
            risk_witness.push_back(r);
            if (a.left.hi <= r && r <= a.right.lo)
                fstar_upper=std::min(fstar_upper,risk(a,pt(r)).hi);
        }
        warm.assign(data.arcs.size(),risk_witness);
    }
    Evaluation evaluate(int h,int k,double q,double r,V eta) {
        ++evaluations;
        const auto& v=data.cache[h][k].variance;
        V x=su(pt(q),data.arcs[h].origin), y=su(pt(r),data.arcs[k].origin);
        V variance=ad(ad(ad(mu(v[0],sq(x)),mu(v[1],mu(x,y))),mu(v[2],sq(y))),
                      ad(ad(mu(v[3],x),mu(v[4],y)),v[5]));
        // Exact variance is nonnegative, even when coefficient cancellation
        // makes its floating interval cross zero.
        variance.lo=std::max(0.0,variance.lo);
        V norm=root(variance);
        Evaluation e;
        e.variance=variance;
        e.value=ad(risk(data.arcs[k],pt(r)),mu(eta,norm));
        if (eta.hi==0) {
            e.gr=scale(su(pt(r),data.arcs[k].mean),2); e.gq=pt(0);
            e.differentiable=true;
        } else if (variance.lo>0 && fc::valid(norm) && std::isfinite(norm.hi)) {
            V vr=ad(ad(scale(mu(v[2],y),2),mu(v[1],x)),v[4]);
            V vq=ad(ad(scale(mu(v[0],x),2),mu(v[1],y)),v[3]);
            e.gr=ad(scale(su(pt(r),data.arcs[k].mean),2),
                    mu(eta,divi(vr,scale(norm,2))));
            e.gq=mu(eta,divi(vq,scale(norm,2)));
            e.differentiable=fc::valid(e.gr)&&fc::valid(e.gq)&&
                std::isfinite(mag(e.gr))&&std::isfinite(mag(e.gq));
        }
        return e;
    }
    // min over r in a closed outer enclosure of the actual chart.
    // Infimum over an interval of possible subgradients is at an endpoint.
    double quadratic_lower(V u,V g) {
        double best=INFINITY;
        for(double z:{g.lo,g.hi}) {
            V half=scale(pt(z),.5);
            best=std::min(best,su(sq(ad(u,half)),sq(half)).lo);
        }
        return best;
    }
    int classify(int h,V q,double cutoff,int iterations=12) {
        const Arc& a=data.arcs[h];
        V eta=divi(pt(cutoff),root(pt(data.n)));
        V diagonal=su(mu(sq(pt(cutoff)),a.variance),scale(sq(su(q,a.mean)),data.n));
        if(diagonal.hi<0) {++diagonal_out;return -1;}
        V fq=risk(a,q);
        if(std::isfinite(fstar_upper)) {
            V fs=pt(fstar_upper);
            double bound=std::min(mu(sq(ad(pt(1),scale(eta,2))),fs).hi,
                ad(fs,mu(scale(pi(),2),mu(eta,root(fs)))).hi);
            if(fq.lo>bound) {++risk_out;return -1;}
        }
        double q0=middle(q);
        V fq0=risk(a,pt(q0)), dq=su(q,pt(q0));
        V lip=ad(pt(fc::up(2*mag(su(q,a.mean)))),scale(mu(eta,root(a.variance)),2));
        double move=mu(lip,pt(mag(dq))).hi;
        bool unknown=diagonal.lo<0;
        bool prepared=false;
        for(size_t k=0;k<data.arcs.size();++k) {
            if(int(k)==h)continue;
            if(risk_lower[k]>=fq.hi) {++risk_skips;continue;}
            if(!prepared){data.prepare(h);prepared=true;}
            ++pairs;
            const Arc& b=data.arcs[k];
            double left=b.left.lo,right=b.right.hi;
            double r=std::max(left,std::min(right,warm[h][k]));
            bool done=false;
            double best_upper_h=INFINITY;
            for(int it=0;it<iterations;++it) {
                Evaluation e=evaluate(h,k,q0,r,eta);
                // A feasible comparison is enough to prove OUT. Lipschitz
                // transports the midpoint inequality across the entire q cell.
                if(b.left.hi<=r && r<=b.right.lo &&
                   su(fq0,e.value).lo>move) {
                    warm[h][k]=r;++rejected;return -1;
                }
                if(e.differentiable) {
                    V u=su(V{b.left.lo,b.right.hi},pt(r));
                    double support_lower=ad(pt(e.value.lo),pt(quadratic_lower(u,e.gr))).lo;
                    double lower=std::max(support_lower,risk_lower[k]);
                    // Joint supporting plane + strong convexity in r:
                    // phi(q,r') >= phi0 + gq*(q-q0) + gr*(r'-r)+(r'-r)^2.
                    double hu=-INFINITY;
                    for(double endpoint:{q.lo,q.hi}) {
                        V support=ad(pt(support_lower),mu(e.gq,su(pt(endpoint),pt(q0))));
                        hu=std::max(hu,su(risk(a,pt(endpoint)),support).hi);
                    }
                    // Independent Lipschitz certificate; either valid bound
                    // suffices. No assumption about convexity of C_n is made.
                    hu=std::min(hu,ad(su(fq0,pt(lower)),pt(move)).hi);
                    best_upper_h=std::min(best_upper_h,hu);
                    if(best_upper_h<=0){warm[h][k]=r;done=true;break;}
                    if(e.gr.lo>0)right=std::min(right,r);
                    else if(e.gr.hi<0)left=std::max(left,r);
                    else {
                        // Strong monotonicity encloses the unique zero.
                        left=std::max(left,su(pt(r),scale(pt(e.gr.hi),.5)).lo);
                        right=std::min(right,su(pt(r),scale(pt(e.gr.lo),.5)).hi);
                    }
                }
                if(!(left<right))break;
                double next=left*.5+right*.5;
                // Newton is a proposal only. The gradient interval above is
                // the only mechanism used to discard parts of the r bracket.
                if(e.differentiable && e.variance.lo>0) {
                    const auto& v=data.cache[h][k].variance;
                    double x=q0-middle(a.origin),y=r-middle(b.origin);
                    double A=middle(v[2]), B=.5*(middle(v[1])*x+middle(v[4]));
                    double C=middle(v[0])*x*x+middle(v[3])*x+middle(v[5]);
                    double vv=middle(e.variance);
                    double curvature=2+middle(eta)*std::max(0.0,A*C-B*B)/std::pow(vv,1.5);
                    double proposal=r-middle(e.gr)/curvature;
                    if(proposal>left && proposal<right && std::isfinite(proposal))next=proposal;
                }
                if(next==r)break;
                r=next;
            }
            warm[h][k]=r;
            if(!done)unknown=true;
        }
        if(unknown){++unresolved;return 0;}
        ++accepted;return 1;
    }
};
}

// Deterministic, full-circle enclosures. Integer cells include every component.
namespace circle_length {
using namespace circle_kernel;
constexpr uint32_t M=1u<<24;
V fraction(uint64_t a,uint64_t b){return a==0?pt(0):divi(pt(double(a)),pt(double(b)));}
using State=std::array<int,3>;
struct Leaf { uint32_t lo,hi; State state; };
struct Context {
    Data data;
    convex_length::Engine convex;
    std::array<double,2> lower,upper;
    double el;
    uint64_t queries=0,visits=0;
    int home=-1;
    Context(NumericVector y,NumericVector l,NumericVector u,double e):data(y),convex(data),lower{l[0],l[1]},upper{u[0],u[1]},el(e){}
    V endpoint(uint32_t i){return scale(pi(),2.0*i/M-1);}
    State classify(const Leaf& leaf) {
        ++queries;
        State s=leaf.state;
        V a=endpoint(leaf.lo),b=endpoint(leaf.hi),q=span(a,b);
        if(home<0 || !(data.arcs[home].left.hi<q.lo && q.hi<data.arcs[home].right.lo)) {
            home=-1;
            for(size_t h=0;h<data.arcs.size();++h)
                if(data.arcs[h].left.hi<q.lo && q.hi<data.arcs[h].right.lo){home=h;break;}
        }
        if(home<0)return s; // True antipodes are enclosed, never assigned a chart arbitrarily.
        for(int j=0;j<2;++j)if(s[j]==0) {
            auto at=[&](double cut){
                int z=convex.classify(home,q,cut);
                return z?z:data.classify(home,q,cut,16,4,visits);
            };
            int l=at(lower[j]),u=l==1?1:(upper[j]==lower[j]?l:at(upper[j]));
            s[j]=l==1?1:u==-1?-1:(l==-1&&u==1)?2:0;
        }
        if(s[2]==0)s[2]=data.classify_el(data.arcs[home],a,b,el);
        // Same membership statistic and nested calibration classes imply C_L subset C_G.
        if(s[1]==1&&s[0]==-1)stop("contradictory G/L certificates");
        if(s[1]==1)s[0]=1;
        if(s[0]==-1)s[1]=-1;
        return s;
    }
};
List measure(NumericVector y,NumericVector lower,NumericVector upper,double el,int budget,double tolerance) {
    fc::check_platform();
    if(y.size()<2||y.size()>200||lower.size()!=2||upper.size()!=2||
       !R_finite(lower[0])||!R_finite(lower[1])||!R_finite(upper[0])||!R_finite(upper[1])||
       lower[0]<lower[1]||upper[0]<upper[1]||lower[1]<0||upper[0]<lower[0]||upper[1]<lower[1]||
       !R_finite(el)||el<0||budget<1||!R_finite(tolerance)||tolerance<=0)
        stop("invalid Circle length controls");
    Context c(y,lower,upper,el);
    struct Priority {bool operator()(const Leaf&a,const Leaf&b)const{return a.hi-a.lo<b.hi-b.lo;}};
    std::priority_queue<Leaf,std::vector<Leaf>,Priority> queue;
    std::vector<uint32_t> cuts{0,M};
    for(const auto&a:c.data.arcs) {
        V index=scale(ad(divi(a.right,pi()),pt(1)),M*.5);
        for(double x:{std::floor(index.lo),std::ceil(index.hi)})
            cuts.push_back(uint32_t(std::max(0.0,std::min(double(M),x))));
    }
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    for(size_t j=1;j<cuts.size();++j)queue.push({cuts[j-1],cuts[j],{0,0,0}});
    std::array<uint64_t,3> unknown{M,M,M},inside{},outside{},band{};
    std::vector<Leaf> leaves;
    int used=0;
    while(!queue.empty()&&used<budget) {
        Leaf a=queue.top();queue.pop();++used;
        State old=a.state;a.state=c.classify(a);
        for(int m=0;m<3;++m)
            if((old[m]==0||old[m]==2)&&(a.state[m]==1||a.state[m]==-1))unknown[m]-=a.hi-a.lo;
        if((a.state[0]==0||a.state[1]==0||a.state[2]==0)&&a.hi-a.lo>1) {
            uint32_t mid=a.lo+(a.hi-a.lo)/2;
            queue.push({a.lo,mid,a.state});queue.push({mid,a.hi,a.state});
        } else leaves.push_back(a);
        if(used%32==0)checkUserInterrupt();
        if(360.0*unknown[0]/M<=tolerance&&360.0*unknown[1]/M<=tolerance&&
           360.0*unknown[2]/M<=tolerance)break;
    }
    while(!queue.empty()){leaves.push_back(queue.top());queue.pop();}
    std::sort(leaves.begin(),leaves.end(),[](const Leaf&a,const Leaf&b){return a.lo<b.lo;});
    IntegerMatrix partition(leaves.size(),5);
    uint64_t total=0;
    for(size_t j=0;j<leaves.size();++j) {
        const auto&a=leaves[j];if(a.lo!=total||a.lo>=a.hi)stop("invalid full-circle partition");total=a.hi;
        partition(j,0)=a.lo;partition(j,1)=a.hi;
        for(int m=0;m<3;++m) {
            partition(j,m+2)=a.state[m];
            if(a.state[m]==1)inside[m]+=a.hi-a.lo;
            if(a.state[m]==-1)outside[m]+=a.hi-a.lo;
            if(a.state[m]==2)band[m]+=a.hi-a.lo;
        }
    }
    if(total!=M)stop("incomplete full-circle partition");
    NumericMatrix bounds(3,2);NumericVector width(3),band_width(3),degrees_width(3);
    for(int m=0;m<3;++m) {
        bounds(m,0)=std::max(0.0,mu(scale(pi(),2),fraction(inside[m],M)).lo);
        bounds(m,1)=mu(scale(pi(),2),fraction(M-outside[m],M)).hi;
        width[m]=su(pt(bounds(m,1)),pt(bounds(m,0))).hi;
        band_width[m]=360.0*band[m]/M;
        degrees_width[m]=360.0*(M-inside[m]-outside[m])/M;
    }
    return List::create(_["bounds"]=bounds,_["width"]=width,_["cutoff_band_width"]=band_width,
        _["degrees_width"]=degrees_width,_["partition"]=partition,_["grid_size"]=double(M),_["inside"]=wrap(inside),
        _["outside"]=wrap(outside),_["band"]=wrap(band),_["queries"]=double(c.queries),
        _["r_cells"]=double(c.visits),_["budget_used"]=used);
}
}

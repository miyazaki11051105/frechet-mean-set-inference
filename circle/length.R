# Circle G/L/EL angular certificates on the run's frozen observations and draws.
# G/L target the full-supremum order statistic of these B draws (not B=infinity).
circle_length_controls <- function(enabled = TRUE, tolerance = .57, stages = 16L,
                                   initial_cells = 8192L, maximum_cells = 65536L,
                                   calibration_cells = 262144L) {
    stopifnot(is.logical(enabled),length(enabled)==1L,!is.na(enabled),
        is.finite(tolerance),tolerance>0,tolerance<=.57,
        stages==as.integer(stages),stages>=0L,stages<=20L,
        initial_cells==as.integer(initial_cells),initial_cells>=256L,
        maximum_cells==as.integer(maximum_cells),maximum_cells>=initial_cells,maximum_cells<=1048576L,
        calibration_cells==as.integer(calibration_cells),calibration_cells>=65536L,calibration_cells<=8388608L)
    list(enabled=enabled,tolerance=tolerance,stages=as.integer(stages),
        initial_cells=as.integer(initial_cells),maximum_cells=as.integer(maximum_cells),
        calibration_cells=as.integer(calibration_cells))
}

circle_length_enclosure <- function(y, lower, upper, el=qchisq(.95,1L),
                                    tolerance=.57, budget=8192L) {
    load_numerics()
    if(length(lower)==1L)lower<-rep(lower,2L)
    if(length(upper)==1L)upper<-rep(upper,2L)
    z <- circle_length_cpp(y,lower,upper,el,as.integer(budget),tolerance)
    p <- z$partition; M <- z$grid_size; d <- p[,2L]-p[,1L]
    stopifnot(M==2^24,p[1L,1L]==0L,tail(p[,2L],1L)==M,all(d>0),
        all(p[-1L,1L]==head(p[,2L],-1L)),all(p[,3:5]%in%c(-1L,0L,1L,2L)),
        !any(p[,5L]==2L),all(is.finite(z$bounds)),all(z$bounds[,1L]>=0),
        all(z$bounds[,1L]<=z$bounds[,2L]))
    for(j in 1:3)stopifnot(z$inside[j]==sum(d[p[,j+2L]==1L]),
        z$outside[j]==sum(d[p[,j+2L]==-1L]),z$band[j]==sum(d[p[,j+2L]==2L]))
    z
}

circle_length_input <- function(output,id,manifest) {
    ip <- file.path(output,'inputs',paste0(id,'.rds'))
    cp <- file.path(output,'certificates',paste0(id,'-circle.rds'))
    input <- readRDS(ip); cert <- readRDS(cp)
    f <- verify_checkpoint(cert,manifest$signature_hash,unname(tools::md5sum(ip)))$fit
    stopifnot(input$family=='circle',f$domain=='full_circle',f$B==1999L,f$k==1900L,
        identical(f$draw_hash,object_hash(f$draws)),
        identical(f$draw_hash,input$original_draw_hash),
        identical(object_hash(input$observations),input$observations_hash),
        identical(f$backend$raw_observations,input$observations))
    list(input=input,fit=f,files=tools::md5sum(c(ip,cp)))
}

circle_length_part <- function(path,signature,files=NULL,paths=NULL) {
    z <- readRDS(path)
    stopifnot(is.character(z$files),length(z$files)==2L,!anyNA(z$files),
        length(names(z$files))==2L,!anyNA(names(z$files)),!anyDuplicated(names(z$files)),
        all(grepl('^[0-9a-f]{32}$',z$files)),
        identical(z$signature,signature),identical(z$payload_hash,object_hash(z$payload)),
        identical(unname(tools::md5sum(names(z$files))),unname(z$files)))
    if(!is.null(files))stopifnot(identical(z$files,files))
    if(!is.null(paths))stopifnot(identical(names(z$files),paths))
    z$payload
}

circle_length_case <- function(output, id, manifest, request) {
    directory <- file.path(output, 'angles', 'cases', id)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    lock <- acquire_lock(file.path(directory, 'case.lock'), 2)
    on.exit(lock$kill(), add = TRUE)
    data <- circle_length_input(output, id, manifest)
    f <- data$fit
    done <- file.path(directory, 'result.rds')
    if (file.exists(done))
        return(circle_length_part(done, request$hash, data$files))
    cfg <- manifest$signature$config$lengths
    started <- elapsed_time()
    checkpoint <- file.path(directory, 'checkpoint.rds')
    history <- list()
    stage0 <- 0L
    if (file.exists(checkpoint)) {
        old <- circle_length_part(checkpoint, request$hash, data$files)
        f$bounds <- intersect_nested_bounds(
            pmax(f$bounds$lower, old$bounds$lower),
            pmin(f$bounds$upper, old$bounds$upper))
        history <- old$history
        stage0 <- old$stage
    }
    stopifnot(stage0 <= cfg$stages)
    f$cutoff <- if (file.exists(checkpoint)) old$cutoff else
        cutoff_order_statistic(f$bounds$lower, f$bounds$upper, f$k, .01)
    el <- unname(qchisq(1 - f$alpha, 1L))
    task <- manifest$signature$plan[
        match(id, manifest$signature$plan$id), , drop = FALSE]
    old_coverage <- evaluate_mean_set(data$fit, task, 'circle')
    write_part <- function(payload, path) save_atomic(
        list(signature = request$hash, files = data$files,
            payload = payload, payload_hash = object_hash(payload)), path)
    for (stage in seq.int(stage0, cfg$stages)) {
        answer <- circle_length_enclosure(
            data$input$observations, f$cutoff$lower, f$cutoff$upper,
            el, cfg$tolerance, cfg$initial_cells)
        if (any(answer$degrees_width > cfg$tolerance) &&
            (answer$degrees_width[3L] > cfg$tolerance ||
                any(answer$cutoff_band_width[1:2] <= cfg$tolerance &
                    answer$degrees_width[1:2] > cfg$tolerance)))
            answer <- circle_length_enclosure(
                data$input$observations, f$cutoff$lower, f$cutoff$upper,
                el, cfg$tolerance, cfg$maximum_cells)
        history_index <- which(vapply(history, function(x) x$stage == stage, logical(1L)))
        stopifnot(length(history_index) <= 1L)
        if (!length(history_index)) history_index <- length(history) + 1L
        history[[history_index]] <- list(
            stage = stage, cutoff = f$cutoff[c('lower', 'upper', 'gap')],
            width = answer$degrees_width, cutoff_band_width = answer$cutoff_band_width,
            queries = answer$queries)
        record <- list(
            id = id, stage = stage, bounds = f$bounds, cutoff = f$cutoff, el = el,
            draw_hash = f$draw_hash, input_hash = data$input$observations_hash,
            history = history, length = answer)
        write_part(record, checkpoint)
        if (all(answer$degrees_width <= cfg$tolerance) || stage == cfg$stages) break
        # Refine only methods whose angular enclosure is still too wide.
        requested <- answer$degrees_width[1:2] > cfg$tolerance
        if (!any(requested)) break
        cells <- as.integer(min(cfg$calibration_cells, max(cfg$initial_cells,
            2 * (length(data$input$observations) + 1L)^2L, 512 * 2^stage)))
        # Avoid expensive coarse Local passes. The complete root cover is already
        # affordable, and the rank solver keeps rigorous bounds on every draw.
        for (j in which(requested)) {
            gap <- f$cutoff$gap[j]
            target <- if (gap < .02) max(.000001, gap / 4) else min(.01, max(.00001, gap / 2))
            policy <- refinement_policy(
                target, stages = 1L, max_cells = cells, max_depth = 13L,
                pair_batch = 64L, draw_batch = 512L,
                direction_tolerance = max(.00001, .05 / 2^stage))
            before <- f$bounds
            f <- refine_cutoff(f, policy, requested_methods = c(j == 1L, j == 2L))
            stopifnot(all(f$bounds$lower >= before$lower),
                all(f$bounds$upper <= before$upper),
                identical(f$draw_hash, object_hash(f$draws)))
        }
        # Save immediately after calibration as well as after measurement.
        record$bounds <- f$bounds
        record$cutoff <- f$cutoff
        record$stage <- stage + 1L
        write_part(record, checkpoint)
    }
    coverage <- evaluate_mean_set(f, task, 'circle')
    stopifnot(all(coverage$coverage_lower >= old_coverage$coverage_lower),
        all(coverage$coverage_upper <= old_coverage$coverage_upper))
    # This audit independently recomputes all 1,999 order statistics with Fraction.
    record$cutoff_audit <- verify_cutoff_gap(f, manifest$signature$config$decimal_gap)
    record$rows <- do.call(rbind, lapply(1:3, function(j) data.frame(
        id = id, n = task$n, parameter = task$parameter, replication = task$replication,
        method = c('G', 'L', 'EL')[j],
        status = if (answer$degrees_width[j] <= cfg$tolerance)
            'ANGLE_PRECISION_ATTAINED' else 'VALID_BOUNDS_BUDGET_LIMIT',
        radian_lower_hex = sprintf('%a', answer$bounds[j, 1L]),
        radian_upper_hex = sprintf('%a', answer$bounds[j, 2L]),
        inside = answer$inside[j], outside = answer$outside[j], band = answer$band[j],
        grid_size = answer$grid_size,
        degrees_lower = 360 * answer$inside[j] / answer$grid_size,
        degrees_upper = 360 * (answer$grid_size - answer$outside[j]) / answer$grid_size,
        degrees_width = 360 * (answer$grid_size - answer$inside[j] - answer$outside[j]) / answer$grid_size,
        degrees_lower_hex = sprintf('%a', 360 * answer$inside[j] / answer$grid_size),
        degrees_upper_hex = sprintf('%a', 360 * (answer$grid_size - answer$outside[j]) / answer$grid_size),
        cutoff_lower = if (j <= 2L) f$cutoff$lower[j] else el,
        cutoff_upper = if (j <= 2L) f$cutoff$upper[j] else el,
        input_md5 = unname(data$files[1L]), draw_hash = f$draw_hash, stage = stage)))
    record$coverage <- data.frame(
        id = id, n = task$n, parameter = task$parameter,
        method = c('G', 'L'), coverage_lower = coverage$coverage_lower,
        coverage_upper = coverage$coverage_upper,
        original_lower = old_coverage$coverage_lower, original_upper = old_coverage$coverage_upper,
        target = 'SAME_DRAW_IDEAL_CUTOFF', cutoff_lower = f$cutoff$lower,
        cutoff_upper = f$cutoff$upper)
    record$elapsed_seconds <- elapsed_time() - started
    stopifnot(identical(unname(tools::md5sum(names(data$files))), unname(data$files)))
    write_part(record, done)
    record
}

circle_length_worker <- function(output) {
    lease<-acquire_lock(file.path(output,'angles','pool.lock'),180,shared=TRUE)
    on.exit(lease$kill(),add=TRUE)
    manifest <- readRDS(file.path(output,'manifest.rds'))
    .libPaths(c(manifest$signature$private_library,.libPaths()))
    .circle_length_worker_main <- function() {
        request <- readRDS(file.path(output,'angles','REQUEST.rds'))
        stopifnot(identical(manifest$signature_hash,object_hash(manifest$signature)),
            identical(unname(tools::md5sum(file.path(output,'source',names(manifest$signature$source_md5)))),
                unname(manifest$signature$source_md5)),identical(request$hash,object_hash(request$signature)))
        load_numerics()
        repeat {
            lock<-acquire_lock(file.path(output,'angles','queue.lock'),60)
            queue<-readRDS(file.path(output,'angles','queue.rds'))
            if(queue$next_index>length(queue$ids)){lock$kill();break}
            id<-queue$ids[queue$next_index];queue$next_index<-queue$next_index+1L
            save_atomic(queue,file.path(output,'angles','queue.rds'));lock$kill()
            circle_length_case(output,id,manifest,request)
        }
    }
    .circle_length_worker_main()
}

circle_length_paths <- function(output)file.path(output,'angles',c(
    'raw.csv','summary.csv','paired.csv','ideal-coverage.csv','table.tex','table.csv','audit.json'))

circle_length_request <- function(manifest) {
    cfg<-manifest$signature$config;plan<-manifest$signature$plan
    stopifnot(identical(manifest$signature_hash,object_hash(manifest$signature)),
        identical(manifest$signature$expected,expected_counts(cfg)),isTRUE(cfg$lengths$enabled),
        identical(cfg$lengths,do.call(circle_length_controls,cfg$lengths)),cfg$B==1999L,cfg$alpha==.05)
    check_report_inventory(plan,simulation_plan(cfg)[c('id','family','parameter','n','replication')],'id')
    signature<-list(schema='circle-global-local-el-angles-v1',display_units='degrees',run=manifest$signature_hash,
        ids=plan$id[plan$family=='circle'],controls=cfg$lengths,
        target='SAME_DRAW_FULL_SUPREMUM_ORDER_STATISTIC',B=cfg$B,k=1900L,EL=sprintf('%a',qchisq(.95,1L)))
    list(signature=signature,hash=object_hash(signature))
}

verify_circle_lengths <- function(output) {
    directory<-file.path(output,'angles');status<-readRDS(file.path(directory,'STATUS.rds'))
    manifest<-readRDS(file.path(output,'manifest.rds'));request<-readRDS(file.path(directory,'REQUEST.rds'))
    expected_request<-circle_length_request(manifest);cfg<-manifest$signature$config
    plan<-manifest$signature$plan;plan<-plan[plan$family=='circle',,drop=FALSE]
    paths<-file.path(directory,'cases',plan$id,'result.rds');reports<-circle_length_paths(output)
    stopifnot(identical(request,expected_request),identical(status$signature,request$hash),
        is.character(status$parts),identical(names(status$parts),paths),
        is.character(status$reports),identical(names(status$reports),reports),
        !anyNA(status$parts),!anyNA(status$reports),
        all(grepl('^[0-9a-f]{32}$',c(status$parts,status$reports))),
        all(c('status','rows','attained','wide')%in%names(status)),
        all(vapply(status[c('rows','attained','wide')],function(x)
            is.numeric(x)&&length(x)==1L&&is.finite(x)&&x>=0&&x==floor(x),logical(1))),
        is.character(status$status),length(status$status)==1L,!is.na(status$status),
        identical(unname(tools::md5sum(names(status$reports))),unname(status$reports)),
        identical(unname(tools::md5sum(names(status$parts))),unname(status$parts)))
    inputs<-file.path(output,'inputs',paste0(plan$id,'.rds'))
    stopifnot(identical(unname(tools::md5sum(inputs)),unname(plan$input_md5)))
    rows<-coverage<-vector('list',nrow(plan))
    for(i in seq_len(nrow(plan))) {
        part<-circle_length_part(paths[i],request$hash,paths=c(inputs[i],
            file.path(output,'certificates',paste0(plan$id[i],'-circle.rds'))))
        stopifnot(identical(part$id,plan$id[i]))
        expected<-cbind(plan[rep(i,3L),c('id','parameter','n','replication')],method=c('G','L','EL'))
        check_report_inventory(part$rows,expected,c('id','method'))
        check_report_inventory(part$coverage,expected[1:2,c('id','parameter','n','method')],c('id','method'))
        stopifnot(setequal(names(part$rows),c('id','n','parameter','replication','method','status',
            'radian_lower_hex','radian_upper_hex','inside','outside','band','grid_size','degrees_lower',
            'degrees_upper','degrees_width','degrees_lower_hex','degrees_upper_hex','cutoff_lower',
            'cutoff_upper','input_md5','draw_hash','stage')),!anyDuplicated(names(part$rows)),!anyNA(part$rows),
            setequal(names(part$coverage),c('id','n','parameter','method','coverage_lower','coverage_upper',
                'original_lower','original_upper','target','cutoff_lower','cutoff_upper')),
            !anyDuplicated(names(part$coverage)),!anyNA(part$coverage),
            all(part$rows$status%in%c('ANGLE_PRECISION_ATTAINED','VALID_BOUNDS_BUDGET_LIMIT')),
            length(part$draw_hash)==1L,is.character(part$draw_hash),grepl('^[0-9a-f]{32}$',part$draw_hash),
            is.numeric(part$stage),length(part$stage)==1L,is.finite(part$stage),
            part$stage==floor(part$stage),part$stage>=0L,part$stage<=cfg$lengths$stages,
            all(part$rows$cutoff_lower>=0),all(part$rows$cutoff_lower<=part$rows$cutoff_upper),
            all(part$coverage$target=='SAME_DRAW_IDEAL_CUTOFF'),
            all(as.matrix(part$coverage[c('coverage_lower','coverage_upper','original_lower','original_upper')])%in%0:1),
            all(part$coverage$original_lower<=part$coverage$coverage_lower),
            all(part$coverage$coverage_lower<=part$coverage$coverage_upper),
            all(part$coverage$coverage_upper<=part$coverage$original_upper),
            all(part$coverage$cutoff_lower>=0),all(part$coverage$cutoff_lower<=part$coverage$cutoff_upper))
        p<-part$length$partition;M<-part$length$grid_size
        stopifnot(is.numeric(M),length(M)==1L,!is.na(M),
            all(vapply(part$length[c('inside','outside','band')],function(x)
                is.numeric(x)&&length(x)==3L&&all(is.finite(x))&&all(x>=0&x==floor(x)),logical(1))),
            is.matrix(p),ncol(p)==5L,nrow(p)>0L,!anyNA(p),all(is.finite(p)),
            all(p==floor(p)),M==2^24,p[1L,1L]==0L,tail(p[,2L],1L)==M,
            all(p[,2L]>p[,1L]),all(p[-1L,1L]==head(p[,2L],-1L)),
            all(p[,3:5]%in%c(-1L,0L,1L,2L)),!any(p[,5L]==2L))
        width<-p[,2L]-p[,1L]
        for(j in 1:3)stopifnot(part$length$inside[j]==sum(width[p[,j+2L]==1L]),
            part$length$outside[j]==sum(width[p[,j+2L]==-1L]),
            part$length$band[j]==sum(width[p[,j+2L]==2L]))
        rows[[i]]<-part$rows;coverage[[i]]<-part$coverage
        index<-match(c('G','L','EL'),part$rows$method)
        stopifnot(all(part$rows$inside[index]==part$length$inside),
            all(part$rows$outside[index]==part$length$outside),all(part$rows$band[index]==part$length$band),
            all(part$rows$grid_size==M),all(part$rows$input_md5==plan$input_md5[i]),
            all(part$rows$draw_hash==part$draw_hash),all(part$rows$stage==part$stage))
    }
    raw<-do.call(rbind,rows);coverage<-do.call(rbind,coverage)
    attained<-sum(raw$status=='ANGLE_PRECISION_ATTAINED');wide<-nrow(raw)-attained
    stopifnot(status$rows==manifest$signature$expected$length_rows,status$rows==nrow(raw),
        status$attained==attained,status$wide==wide,
        status$status==if(wide==0L)'COMPLETE_ANGLE_PRECISION_ATTAINED'else'COMPLETE_WITH_VALID_WIDE_ANGLE_BOUNDS')
    # Rebuild reports only in a temporary directory. The existing exact-arithmetic
    # auditor checks the published integer measures, precision and rational means.
    temporary<-tempfile('angle-audit-');dir.create(temporary)
    on.exit(unlink(temporary,recursive=TRUE),add=TRUE)
    save_atomic(raw,file.path(temporary,'raw.csv'));save_atomic(coverage,file.path(temporary,'ideal-coverage.csv'))
    keys<-list(raw=c('id','method'),`ideal-coverage`=c('id','method'),
        summary=c('parameter','n','method'),paired=c('parameter','n','contrast'))
    compare_csv<-function(name)stopifnot(identical(
        reproduction_csv_values(file.path(directory,paste0(name,'.csv')),keys[[name]]),
        reproduction_csv_values(file.path(temporary,paste0(name,'.csv')),keys[[name]])))
    compare_csv('raw');compare_csv('ideal-coverage')
    program<-normalizePath(file.path('common','exact.py'),mustWork=TRUE)
    if(system2('python3',c(shQuote(program),'angles',shQuote(temporary),as.character(cfg$reps),
        sprintf('%a',cfg$lengths$tolerance)))!=0L)stop('Independent read-only angular audit failed')
    compare_csv('summary');compare_csv('paired')
    stopifnot(all(file.copy(file.path(directory,c('table.csv','table.tex')),temporary)))
    if(system2('python3',c(shQuote(program),'angle-table',shQuote(temporary)))!=0L)
        stop('Independent read-only angular table audit failed')
    stopifnot(identical(jsonlite::read_json(file.path(directory,'audit.json')),
        jsonlite::read_json(file.path(temporary,'audit.json'))))
    invisible(status)
}

write_circle_length_table <- function(summary,directory,pilot) {
    rows<-list();lines<-c('% All angles are in degrees; outward decimal endpoints.',
        '% Numerical enclosures of realized averages, not Monte Carlo confidence intervals.',
        '\\begin{tabular}{rrccc}',
        sprintf('\\multicolumn{5}{c}{%s$R=%d$; degrees} \\\\',if(pilot)'PILOT; 'else'',summary$replications[1]), '$\\lambda$ & $n$ & Global & Local & EL \\\\', '\\hline')
    keys<-unique(summary[c('parameter','n')]);keys<-keys[order(keys$parameter,keys$n),]
    for(i in seq_len(nrow(keys))) {
        z<-summary[summary$parameter==keys$parameter[i]&summary$n==keys$n[i],]
        text<-vapply(c('G','L','EL'),function(m){a<-z[z$method==m,]
            paste0('[',a$mean_degrees_lower_tex,', ',a$mean_degrees_upper_tex,']')},character(1L))
        rows[[i]]<-data.frame(lambda=keys$parameter[i],n=keys$n[i],Global=text[1],Local=text[2],EL=text[3],replications=z$replications[1])
        lines<-c(lines,sprintf('%g & %d & %s & %s & %s \\\\',keys$parameter[i],keys$n[i],text[1],text[2],text[3]))
    }
    writeLines(c(lines,'\\end{tabular}'),file.path(directory,'table.tex'))
    table<-do.call(rbind,rows);save_atomic(table,file.path(directory,'table.csv'))
}

run_circle_lengths <- function(output, project) {
    manifest <- readRDS(file.path(output, 'manifest.rds'))
    cfg <- manifest$signature$config
    if (!isTRUE(cfg$lengths$enabled)) return(invisible(NULL))

    directory <- file.path(output, 'angles')
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    # A worker keeps a shared lease even if its controller disappears. A new
    # controller must obtain the exclusive guard before rewriting the queue.
    guard <- acquire_lock(file.path(directory, 'pool.lock'), 2)
    on.exit({if (!is.null(guard)) guard$kill()}, add = TRUE)
    plan <- manifest$signature$plan
    ids <- plan$id[plan$family == 'circle']
    request <- circle_length_request(manifest)
    request_path <- file.path(directory, 'REQUEST.rds')
    if (file.exists(request_path)) stopifnot(identical(request, readRDS(request_path))) else
        save_atomic(request, request_path)
    status_path <- file.path(directory, 'STATUS.rds')
    if (file.exists(status_path)) {
        previous <- readRDS(status_path)
        stopifnot(identical(previous$signature, request$hash))
        present <- file.exists(names(previous$reports))
        stopifnot(identical(unname(tools::md5sum(names(previous$reports)[present])),
            unname(previous$reports[present])))
        if (all(file.exists(circle_length_paths(output)))) return(verify_circle_lengths(output))
    }
    paths <- file.path(directory, 'cases', ids, 'result.rds')
    existing <- file.exists(paths)
    for (path in paths[existing]) circle_length_part(path, request$hash)
    save_atomic(list(ids = ids[!existing], next_index = 1L), file.path(directory, 'queue.rds'))
    workers <- min(cfg$workers, sum(!existing))
    processes <- list()
    on.exit({for (p in processes) if (p$is_alive()) p$kill_tree()}, add = TRUE)
    if (workers > 0L) for (w in seq_len(workers)) {
        processes[[w]] <- processx::process$new(
            file.path(R.home('bin'), 'Rscript'),
            c('--vanilla', file.path(project, 'run.R'), 'length-worker', output), wd = project,
            stdout = file.path(directory, paste0('worker-', w, '.log')),
            stderr = file.path(directory, paste0('worker-', w, '.err')), cleanup_tree = TRUE)
    }
    guard$kill()
    guard <- NULL
    repeat {
        if (any(vapply(processes, function(p)
            !p$is_alive() && !identical(p$get_exit_status(), 0L), logical(1L))))
            stop('Circle angle worker failed; inspect angles/worker-*.err and rerun to resume verified checkpoints.')
        if (!length(processes) || !any(vapply(processes, function(p) p$is_alive(), logical(1L)))) break
        Sys.sleep(10)
    }
    failed <- which(vapply(processes, function(p) p$get_exit_status() != 0L, logical(1L)))
    if (length(failed))
        stop('Circle angle worker failed; checkpoints retained. Inspect angles/worker-*.err; rerun to resume.')
    stopifnot(all(file.exists(paths)))
    parts <- lapply(paths, circle_length_part, signature = request$hash)
    raw <- do.call(rbind, lapply(parts, `[[`, 'rows'))
    coverage <- do.call(rbind, lapply(parts, `[[`, 'coverage'))
    stopifnot(nrow(raw) == 36L * cfg$reps, !anyDuplicated(paste(raw$id, raw$method)))
    save_atomic(raw, file.path(directory, 'raw.csv'))
    save_atomic(coverage, file.path(directory, 'ideal-coverage.csv'))
    result <- system2('python3', c('common/exact.py', 'angles', shQuote(directory),
        as.character(cfg$reps), sprintf('%a', cfg$lengths$tolerance)))
    if (result != 0L) stop('Independent Circle angle/report arithmetic audit failed')
    summary <- read.csv(file.path(directory, 'summary.csv'),
        colClasses = c(mean_degrees_lower_tex = 'character', mean_degrees_upper_tex = 'character'))
    write_circle_length_table(summary, directory, cfg$reps != 1000L)
    if (system2('python3', c('common/exact.py', 'angle-table', shQuote(directory))) != 0L)
        stop('Independent angular table audit failed')
    status <- list(
        status = if (all(raw$status == 'ANGLE_PRECISION_ATTAINED'))
            'COMPLETE_ANGLE_PRECISION_ATTAINED' else
            'COMPLETE_WITH_VALID_WIDE_ANGLE_BOUNDS',
        signature = request$hash, rows = nrow(raw),
        attained = sum(raw$status == 'ANGLE_PRECISION_ATTAINED'),
        wide = sum(raw$status != 'ANGLE_PRECISION_ATTAINED'),
        parts = tools::md5sum(paths), reports = tools::md5sum(circle_length_paths(output)))
    save_atomic(status, status_path)
    verify_circle_lengths(output)
    invisible(status)
}

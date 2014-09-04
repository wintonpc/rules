;; utils

(define-syntax (define-memoized stx)
  (define (key-syntax stx)
    (syntax-case stx ()
      [(arg0) #'arg0]
      [(arg0 arg1 arg2 ...) #'(list arg0 arg1 arg2 ...)]))
  (syntax-case stx ()
    [(_ (name . args) exp0 exp1 ...)
     #`(define name
         (let ([memos (make-hash)])
           (lambda args
             (hash-ref! memos #,(key-syntax #'args)
                        (lambda () exp0 exp1 ...)))))]))

(define (pipe input . procs)
  (foldl (lambda (p acc) (p acc)) input procs))

(define-syntax (pipe-each stx)
  (syntax-case stx (in)
    [(_ var in xs tx ...)
     #'(map (lambda (var) (pipe var tx ...)) xs)]))

(define-syntax (foreach stx)
  (syntax-case stx (in)
    [(_ var in xs exp0 exp1 ...)
     #'(for-each (lambda (var) exp0 exp1 ...) xs)]))

(define-syntax (map-each stx)
  (syntax-case stx (in)
    [(_ var in xs exp0 exp1 ...)
     #'(map (lambda (var) exp0 exp1 ...) xs)]))

(define (~> proc . early-args)
  (lambda late-args
    (apply proc (append early-args late-args))))

(define first car)

(define (non-false x) x)

(define cons*
  (case-lambda
    [(xs) xs]
    [(h . t) (cons h (apply cons* t))]))

(define-syntax (is-match? stx)
  (syntax-case stx ()
    [(_ x pat ...)
     #'(match x
         [pat #t]
         ...
         [_ #f])]))


;; model

(struct Assay (name comp-meths))
(struct CompoundMethod (name settings enabled-for-rule? chrom-meths assay))
(struct ChromatogramMethod (id classifier settings compound-method))

(struct Batch (name assay samples))
(struct Sample (id name type injection dilution compounds [batch #:mutable]))
(struct Compound (id name nom-conc std-dev deviation response use-record? calc-conc
                     attempted-regression? regression-succeeded? found-concentration?
                     chromatograms [sample #:mutable]))
(struct Chromatogram (classifier peak-area expected-rt [compound #:mutable]))

(struct RuleInfo (id description) #:transparent)

(struct Result (rule-id name handle value) #:transparent)
(struct Handle (sample compound chromatogram) #:transparent)
(struct Flag (id text abbreviation category severity) #:transparent)
(struct Flagged (handle flag rule-id) #:transparent)

;; model accessors

(define-memoized (get-comp-meths batch)
  (pipe batch
    (~> Batch-assay)
    (~> Assay-comp-meths)))

(define-memoized (get-comps-by-name batch comp-meth-name)
  (pipe batch
    Batch-samples
    (~> map (~> get-sample-comp comp-meth-name))))

(define (is-compound-named comp-name)
  (lambda (comp)
    (string=? (Compound-name comp) comp-name)))

(define-memoized (get-sample-comp comp-name sample)
  (pipe sample
    (~> Sample-compounds)
    (~> findf (is-compound-named comp-name))))

(define (standard-compound? comp)
  (standard-sample? (Compound-sample comp)))

(define-memoized (standard-sample? samp)
  (string=? (Sample-type samp) "standard"))

(define (get-quant-chrom comp)
  (let ([chroms (filter is-quant-chrom (Compound-chromatograms comp))])
    (if (null? chroms) #f chroms)))

(define (is-quant-chrom chrom)
  (string=? (Chromatogram-classifier chrom) "Quant")) ; fix

(define (mean xs)
  (if (null? xs)
      +nan.0
      (/ (foldl + 0 xs) (length xs))))

(define-memoized (get-compound-names batch)
  (map CompoundMethod-name (get-comp-meths batch)))



;; rules engine

(define (run-rules batch . rules)
  (define (run flags results rules)
    (if (null? rules)
        (values flags results)
        (let* ([rule (car rules)]
               [info (rule 'reflect)]
               [rule-id (RuleInfo-id info)])
          (let-values ([(fs rs) (rule 'run batch
                                      (lambda (n h v) (Result rule-id n h v))
                                      (~> lookup-result results))])
            (run (append flags fs) (append results rs) (cdr rules))))))
  (run '() '() rules))


(define (lookup-result results result-name compound-name)
  (first
   (filter 
    (lambda (r)
      (is-match? r (Result _ (== result-name) (Handle #f (== compound-name) #f) _)))
    results)))

(define (reflect-rules . rules)
  (map (lambda (rule) (rule 'reflect)) rules))

(define (load-batch)
  (load "/home/pwinton/git/rules/batch.rkt")
  (connect-parents the-batch))

(define (connect-parents batch)
  (foreach samp in (Batch-samples batch)
    (set-Sample-batch! samp batch)
    (foreach comp in (Sample-compounds samp)
      (set-Compound-sample! comp samp)
      (foreach chrom in (Compound-chromatograms comp)
        (set-Chromatogram-compound! chrom comp)))))

(define (make-compound-handle comp-name)
  (Handle #f comp-name #f))

;; rules

(define-syntax (define-rule stx)
  (syntax-case stx ()
    [(_ (name batch make-result read-result) desc exp0 exp1 ...)
     #'(define (name . L)
         (match L
           ['(reflect) (RuleInfo 'name desc)]
           [`(run ,batch ,make-result ,read-result)
            (let ()
              exp0 exp1 ...)]))]))

(define-rule (mean-peak-area batch make-result read-result)
  "Calculates the mean quant peak area of standards."
  (let ([results
         (pipe-each comp-name in (get-compound-names batch)
           (~> get-comps-by-name batch)
           (~> filter standard-compound?)
           (~> map get-quant-chrom)
           (~> filter non-false)
           (~> map Chromatogram-peak-area)
           mean
           (lambda (r) (make-result 'mean-peak-area (make-compound-handle comp-name) r)))])
    (values '() results)))

(define-rule (peak-area-deviation batch make-result read-result)
  "blah blah blah"
  (map-each comp-meth in (get-comp-meths batch)
    (let ([mean-area (read-result 'mean-peak-area)]
          [threshold (read-parameter comp-meth "ISAreaDeviationThreshold")]
          [quant-chrom (get-quant-chrom)])
      (let ([deviation ()])
          
      
      
        (values flags '())))))
  
(pretty-print-columns 120)
(load-batch)

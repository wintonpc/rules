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

(define (~> proc . early-args)
  (lambda late-args
    (apply proc (append early-args late-args))))

(define first car)


;; model

(struct <assay> (name comp-meths))
(struct <compound-method> (name settings enabled-for-rule? chrom-meths assay))
(struct <chromatogram-method> (id classifier settings compound-method))

(struct <batch> (name assay samples))
(struct <sample> (id name type injection dilution compounds [batch #:mutable]))
(struct <compound> (id name nom-conc std-dev deviation response use-record? calc-conc
                     attempted-regression? regression-succeeded? found-concentration?
                     chromatograms [sample #:mutable]))
(struct <chromatogram> (classifier peak-area expected-rt [compound #:mutable]))


;; model accessors

(define-memoized (get-comp-meths batch)
  (pipe batch
    (~> <batch>-assay)
    (~> <assay>-comp-meths)))

(define-memoized (get-comps-by-name batch comp-meth-name)
  (pipe batch
    <batch>-samples
    (~> map (~> get-sample-comp comp-meth-name))))

(define (is-compound-named comp-name)
  (lambda (comp)
    (string=? (<compound>-name comp) comp-name)))

(define-memoized (get-sample-comp comp-name sample)
  (pipe sample
    (~> <sample>-compounds)
    (~> findf (is-compound-named comp-name))))

(define (standard-compound? comp)
  (standard-sample? (<compound>-sample comp)))

(define-memoized (standard-sample? samp)
  (string=? (<sample>-type samp) "standard"))

(define (get-quant-chroms comp)
  (filter is-quant-chrom (<compound>-chromatograms comp)))

(define (is-quant-chrom chrom)
  (string=? (<chromatogram>-classifier chrom) "Quant")) ; fix

(define (mean xs)
  (if (null? xs)
      +nan.0
      (/ (foldl + 0 xs) (length xs))))

(define-memoized (get-compound-names batch)
  (map <compound-method>-name (get-comp-meths batch)))



;; rules engine

(define (run-rules batch . rules)
  (define (run flags results rules)
    (if (null? rules)
        (values flags results)
        (let-values ([(name fs rs) ((car rules) batch #f)])
          (let ([annotated-rules (map (~> list name) rs)])
            (run (append flags fs) (append results annotated-rules) (cdr rules))))))
  (run '() '() rules))

(define (load-batch)
  (load "/home/pwinton/git/rules/batch.rkt")
  (connect-parents the-batch))

(define (connect-parents batch)
  (foreach samp in (<batch>-samples batch)
    (set-<sample>-batch! samp batch)
    (foreach comp in (<sample>-compounds samp)
      (set-<compound>-sample! comp samp)
      (foreach chrom in (<compound>-chromatograms comp)
        (set-<chromatogram>-compound! chrom comp)))))

;; rules

(define (mean-peak-area batch read-result)
  (let ([results
         (pipe-each comp-name in (get-compound-names batch)
           (~> get-comps-by-name batch)
           (~> filter standard-compound?)
           (~> map get-quant-chroms)
           flatten
           (~> map <chromatogram>-peak-area)
           mean
           (~> cons comp-name))])
    (values 'mean-peak-area '() results)))




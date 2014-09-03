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

(define (|| input . procs)
  (foldl (lambda (p acc) (p acc)) input procs))

(define (~> proc . early-args)
  (lambda late-args
    (apply proc (append early-args late-args))))

(define first car)


;; model

(struct <assay> (name comp-meths))
(struct <compound-method> (name settings enabled-for-rule? chrom-meths))
(struct <chromatogram-method> (id classifier settings))

(struct <batch> (name assay samples))
(struct <sample> (id name type injection dilution compounds))
(struct <compound> (id name nom-conc std-dev deviation response use-record? calc-conc
                     attempted-regression? regression-succeeded? found-concentration?
                     chromatograms))
(struct <chromatogram> (classifier peak-area expected-rt))


;; model accessors

(define-memoized (get-comp-meths batch)
  (|| batch (~> <batch>-assay) (~> <assay>-comp-meths)))

(define-memoized (get-comp batch comp-meth-name)
  (|| batch <batch>-samples (~> map (~> get-sample-comp comp-meth-name))))

(define-memoized (get-sample-comp comp-name sample)
  (|| sample (~> <sample>-compounds) (~> findf (lambda (comp) (string=? (<compound>-name comp) (comp-name))))))


;; rules

(define (mean-peak-area batch read-result)
  (|| (get-comp-meths batch)
   (~> map comp-meth-name)
   (~> map (~> get-comp batch))))
            


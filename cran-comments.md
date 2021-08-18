## Test environments

* Windows 10 x64, R 4.0.4
* Mac x64, R 4.0.3

## R CMD check results

  Note: information on .o files for x64 is not available
  File 'C:/Users/damon/Desktop/fMRI/fMRIscrub.Rcheck/fMRIscrub/libs/x64/fMRIscrub.dll':
    Found 'abort', possibly from 'abort' (C), 'runtime' (Fortran)
    Found 'exit', possibly from 'exit' (C), 'stop' (Fortran)
    Found 'printf', possibly from 'printf' (C)
  
  Compiled code should not call entry points which might terminate R nor
  write to stdout/stderr instead of to the console, nor use Fortran I/O
  nor system RNGs. The detected symbols are linked into the code but
  might come from libraries and not actually be called.
  
  See 'Writing portable packages' in the 'Writing R Extensions' manual.
   WARNING
  'qpdf' is needed for checks on size reduction of PDFs

0 errors v | 1 warning x | 1 note x

`abort`, `exit` and `printf` is never called by our code.

## Downstream dependencies

None.

## Tests

Passes all the tests in `tests/testthat.R`

## Resubmission

Please always explain all acronyms in the description text.

* Done!

If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

* Added DOI and arXiv links in this format!

Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation. Please write about the
structure of the output (class) and also what the output means. (If a
function does not return a value, please document that too, e.g.
\value{No return value, called for side effects} or similar)
Missing Rd-tags:
      data_CompCor_params.Rd: \value
      fsl_bptf.Rd: \value
      noise_Params.Rd: \value
      pscrub_Params.Rd: \value
      rob_scale.Rd: \value
      summary.scrub_DVARS.Rd: \value
      summary.scrub_FD.Rd: \value
      summary.scrub_projection.Rd: \value

* Done for all exported methods! (`*_Params` are internal)

You have examples for unexported functions.
Please either omit these examples or export these functions.
Used ::: in documentation:
      man/pscrub_multi.Rd:
         psx = fMRIscrub:::pscrub_multi(X)

* Removed example for unexpored function!

Please always make sure to reset to user's options(), working directory
or par() after you changed it in examples and vignettes and demos.
e.g.:
oldpar <- par(mfrow = c(1,2))
...
par(oldpar)
e.g. inst/doc/projection_scrubbing.R

* Fixed in the vignette!

Please do not set a seed to a specific number within a function.

* We added an argument `seed` to each applicable function to control whether
  the seed is set or not, and if so, what value it is set to.

## Second resubmission

   Possibly misspelled words in DESCRIPTION:
     incudes (27:38)

* Fixed!
# ADiGator Version 1.4

A Source Transformation via Operator Overloading Tool for the Automatic Differentiation of Mathematical Functions Defined by MATLAB Code

Please see doc/ADiGatorUserGuide.pdf for User's Guide.

Copyright 2011-2017 Matthew J. Weinstein and Anil V. Rao
Distributed under the GNU General Public License version 3.0
Please see COPYING.txt for full License.

## Contact Info:
weinstein87@gmail.com

# Citing ADiGator:
Please cite the most recent ACM-TOMS CALGO article. BibTex is here:
```
@article{weinstein2017algorithm,
  title={Algorithm 984: ADiGator, a toolbox for the algorithmic differentiation of mathematical functions in MATLAB using source transformation via operator overloading},
  author={Weinstein, Matthew J and Rao, Anil V},
  journal={ACM Transactions on Mathematical Software (TOMS)},
  volume={44},
  number={2},
  pages={21},
  year={2017},
  publisher={ACM}
}
```

# Release Notes:
V 1.5 6/02/19
* A bit of repo re-organization, moved everything to github: https://github.com/matt-weinstein/adigator
* Fixed derivative rules for a few unary trig functions. Many "degree" versions of trig routines had a bug, and a few of the hyperbolic versions were also off.

V 1.4 4/22/17
* Fix for matlab 2016b changing matlab.codetools.requiredFilesAndProducts to also return .mat files. Was causing issues with 2nd derivatives.
* General code cleanup of comments and examples.
* Added adigatorInput class rather than using cada objects to define inputs - only used to define derivative/auxiliary adigator input variables.
* Added adigatorFindMatchingParen.m - finds matching pair for a parenthesis, curly or square bracket. Modified some old code to use this for parsing function I/O and also for finding multiple commands written on the same line. Users likely experienced some issues when performing multiple catonations on the same line (either erroring out while attempting to generate a temporary function, or getting an invalid expression statement within the temporary function), particularly if doing so on an input to a function. This should fix those issues.
* Fixed a couple of minor bugs for cada/cross and cada/transpose in the vectorized mode.
* Made some changes to handle dead code better, particularly loops which never run, and sub-functions which never get called.
* Changed check in adigatorPrintTempFiles to print out cada#input# and cada#output# assignments in temp functions based upon names read from user file rather than derivative number - was a bug where sub-function numbers change if mixing different adigator generated files, should fix that.

V 1.3 8/21/16
* Fixed a couple of bugs using vectorized repmat.
* Fixed a bug where horzcat/vertcat were not handling numeric inputs properly when called from within a loop.
* Added checks to horzcat/vertcat to just return input if they get only a single input. I found some odd behavior where Matlab seemed to be forcing this behavior.
* Added in some functionality for complex numbers. If you want to use it, there is an options.complex = 1 that needs to be set. It defaults to 0 so that abs/transpose can be simplified (assuming no complex numbers exist). I didn't play with this a ton, but did some baseline validation of complex functions of real inputs.
* Put in a check to ignore "called functions" which have been made into @cada methods. If you have some utility functions which get called repeatedly within a program, it may be more efficient to "overload" them, similar to how I have done dot(), fliplr(), others. Note: if you do this, you need to be careful if the function calls size(), length(), or something else which has identically zero derivatives - there is a workaround in the overloaded fliplr().

V 1.2 2/29/16
* This release is mainly just a code reorganization and cleanup.
* I switched over to Matlab's newer "classdef" class definitions. This allowed me to get rid of a bunch of .m files within the @cada directory and place the routines directly within @cada/cada.m.
* The folder structures have changed a bit, so users will need to run startupadigator.m to set new paths. This was to separate out the four main functions, adigator, adigatorOptions, adigatorCreateDerivInput, adigatorCreateAuxInput, where these used to be contained within the util directory. There is also now a cadaUtils directory within lib/ that contains some routines I used to have stored in multiple private directories.
* I did a pretty thorough sweep through the code and removed some old files/code that wasn't being used anymore.
* The one bug that I fixed was only showing up when performing three or more differentations, with different variables of differentation. For instance, if computing d(d(df/dx)/dx)/dy, the variable names were getting messed up in the f_xxy file. This should be fixed.

V 1.1.1 11/2/15 mjw
* Removed a "keyboard" in @cada/cadaPrintReMap.m
* The newest version of Matlab (2015b) has disallowed the use of depfun, which adigator was using to determine a user called functions list. It is now forcing one to use"matlab.codetools.requiredFilesAndProducts"
  * I changed the code to use this new function if using 2015b or later.
Unfortunately, the new function is extremely slow (roughly 2seconds for the call) which is quite annoying, particularly when differentiating smaller functions which should take only a fraction of a second to do the transformation.
  * If anyone has any insight to a better work-around, or has found a way to use the old depfun, please let me know.

V 1.1 5/8/15 mjw
* Made some changes to the way that inputs/outputs of external called functions are handled. The primary change was to create different input/output naming schemes for each function. So, now if you call 2 functions, myfun1, myfun2, each with a single input, the inputs will now be input1_1, input1_2, whereas before they were both called input1. This was causing some issues when handling auxiliary data and should now be fixed.
* Modified handling of while loops. When I originally coded while loops, I did them wrong and for some reason assumed three iterations would be enough to find a static (in regards to sparsity) loop iteration. These were re-worked to (during code generation) analyze the loops recursively until a static iteration is reached. There is a maximum iteration limit which the user can set.
* Added some rudimentary compression files, adigatorColor and adigatorUncompressJac. If, for some reason, you know derivative sparsity beforehand, utilizing compression can help to reduce file generation times (particular for large variable problems which perform matrix multiplications/summations).
* Modified adigatorGenFiles4Ipopt and adigatorGenFiles4Fmincon to (1) allow for linear objective functions and (2) utilize compression when computing d^2(lambda.'c)/dx^2, where c is a vector of constraints. That is, the columns of the m by nn matrix d^2c/dx^2 are removed when pre-multiplying by lambda.'.
* Minor bug fixes for cadastruct and cada subsasgn and subsref when the references/assignments are loop iteration-dependent or external function call-dependent.
* mtimes, when called within a loop, will now check to see if the summation dimension changes and remove false values if so.
* cleaned up the directories, have all of the adigator transformation routines placed in the lib folder, any functions that the user might call placed in the util folder. Also added Contents files for each folder.
edit 5/20/15:
* There were a couple of bugs in adigatorForIterEnd>CollectChildData when collecting nested loop data of structure/cell references and assignments. Also found a bug in there that would only surface if there is an organizational operation within an iteration dependent conditional statement within a nested loop which runs for differing number of iterations. I fixed what I could find, updated the version number in a couple of places and just replaced the V1.1 release.
edit 6/18/15:
* There was a bug when summing over the second dimension of a matrix, ex: sum(x,2). It should be fixed and updated release.

V 1.0 12/4/14 mjw
Big changes:
* Cell/structure handling - cells and structures are now treated as overloaded objects with their own class (@cadastruct). Previously, only the elements of cell/structure arrays which corresponded to numeric objects were being tracked - this lead to some issues when cell/structure arrays were used with flow control (really an issue for GPOPS2 users for multi-phase problems). Now, if a structure is a scalar, each field is treated as a variable, otherwise all cell/structures are treated as variables. This should take care of any issues people were having when having issues with multiphase problems. It is also a much more elegant solution than the previous.
* Addition of "black box" file generation commands: I realized it was a bit of a pain to always be writing wrapper files for simple problems so I created some routines which both call the adigator command, and create wrapper files for Jacobians, Hessians, gradients, Lagrangrian Hessians, etc. These new files are as follows: adigatorGenJacFile, adigatorGenHesFile, adigatorGenFiles4Fsolve, adigatorGenFiles4Fminunc, adigatorGenFiles4Fmincon, adigatorGenFiles4Ipopt. See the updated user's guide or /help for more information.
* Tweaking of derivative routine for matrix operations: When printing derivatives of matrix operations (e.g. mtimes, sum, etc.), I project into a derivative matrix (rather than operating on vectors of non-zeros as is done for other operations). If a matrix A(x), was [m p] size and x of length n, then I was always projecting into a derivative matrix of size [m p*n] and operating on it from the left. I noticed that for many problems, however, that many columns of the m by pn matrix were known to be zero, particularly at the second derivative level. I now determine which q<=pn columns have non-zero elements, and project into a matrix of size [m q]. This is sort of a very quick coloring scheme (looking into a better coloring scheme might be a good idea and would allow you to take into account the sparsity pattern of function matrices). Moreover, for the b = sum(A(x)) function, I do a check to see if nnz(dA/dx) = q, if so, then no projection is made as there is only single non-zero derivative in each column of the [m q] matrix (think of summing the identity matrix).
* While loops: Certain while loops are now allowed, specifically those which do not contain any iteration dependent operations. This is mainly to allow for performing Newton-type iterations. There is more info in the updated user’s guide.
* User’s guide & examples: I added a couple more examples and modified some others to use the new file generation commands. There is more information on what examples use what adigator commands in the updated user’s guide.
New Overloaded Functions:
* numel - this was requested, it should work properly, but be careful putting keyboards in it
* nan - as requested
* prod - prod was a bit tricky, unfortunately I only have it coded to take the prod of a vector if the input has derivatives, the matrix version adds a dimension which makes it too complicated to concisely print derivatives to file. Even the code that is printed by the vector version is quite cumbersome. That being said, one can now do the matrix prod with one loop instead of two!
* nnz - not sure why I didn’t have this - added it as I used it for writing prod derivs.
* mrdivide - this wasn't set up for the non-square case, I was lazy and just made it call mldivide with some transposes
Bug fixes:
* transpose was not computing function sparsity properly
* horzcat/vertcat were not handling empty inputs properly
* fixed some issues with variables not being pre-allocated prior to loops - NOTE: if you are using vectorized mode (GPOPS2 continuous functions), you should ALWAYS pre-allocate vectorized variables, if you don’t, adigator doesn’t know the vectorized dimension when pre-allocating and this can cause some errors.
* were some errors being generated with adigatorPrintTempFiles was calling genErrorLink routine, these should all be fixed
* made it so that you can have 
```
y = calcs ... %comment
morecalcs; %without error
```
* Update 12/10/14 -- mtimes derivs had a minor bug, wasn't projecting properly when result of c(x)=a(x)*b(x) was a scalar and a,b vectors.

V0.4.4 9/12/14 mjw
* Fixed an issue with logical assignment of scalar to an array, e.g., y(logical) = 0; - was producing errors at second derivative of functions with sqrt in them due to the added derivative check of V0.4.3.
* Changed deriv check of sqrt(x) to check that both x and dx are zero, rather than just x.
* Added similar check for power function x^a, when a < 1.
* I believe there are still some issues pertaining to cell/structure array handling, particularly when referencing/assigning to structure array elements within a loop. The current way in which these are handled is a bit of a kluge and I am currently working on overloading cells and structures so that I can handle this in a better fashion.
EDIT 11/16/14
* Apologies, I apparently got one of the bugs associated with y(logical) = 0, but missed another. I updated the release on here.

V0.4.3 8/13/14 mjw
* Overloaded the following: min, max, fix, floor, round, ceil, rem, mod
* Added a check to the sqrt function such that y = sqrt(x) evaluated at dx = 0, x = 0, produces dy = 0
* Changed up logical referencing/assignment a bit, it was not working properly when called from within loops and it should now
* Fixed an issue with cell/structure array references within loops where an operation had to be performed to get the reference index (e.g. y = s(i+1).y)

V0.4 6/12/14 mjw
* Projecting using the sparse(i,j,x,m,n) command, when x had known zero function locations was messing up, fixed this
* Was an error in a subsref check causing it to yell at you about switching between logical and linear indexing within a loop when you actually weren’t doing that, fixed this
* Was something with re-mapping operation being printed in conditional branches which never fire, this should be fixed
There were a few more changes with the ever difficult structure/cell array referencing/assignments within loops:
  * If had something like si = struct(i).a; and si had a field assigned a string then stuff was messing up, this should be fixed, same with assignment version
  * If you have a statement like
```
for i =1:n
    si = s(i);
end
```
where s is a structure array, I can't allow for you to keep the loop rolled in the derivative file, since the reference looks like s is an array and I have to look for cell/structure references from the source code level. This needs to be replaced by "si = s(i).a" or "si = s{i}" (need to change the way you store stuff), if you attempt to do it the "wrong" way an error message should be produced which explains this. Same goes for assignment counterpart.
  * If you had
```
s = struct('f',cell(3,1));
for i = 1:3
    s(i).f = somestructure;
end
```
and pre-allocated like a good person, the assignment was getting messed up, fixed this for all test cases I tried.
* Regarding interp1 and ppval - I made it so that you can do something like
```
for i = 1:n
    y(i).val = ppval(data(i).pp,x);
end
```
where the pp changes on each iteration. However, the .dim and .order of each pp must be the same. Additionally, if the number of pieces changes, this will result in some redundant calculations being printed to file at the 2nd and higher derivative level (this is kind of a kluge). If it is the case that you have
```for i = 1:n
    y(i).val = interp1(data(i).x,data(i).y,x);
end
```
it will not let you keep this rolled, will need to apriori generate the polys and make them inputs to your function.
* Regarding interp2 - the way I handle interp2 is that I actually generate the 2D piecewise polynomials and write the derivative code in terms of the adigatorEvalInterp2pp command. This command is equivalent to using ppval except MATLAB has no way of generating the coefficients automatically. In this release I allowed for you to create a 2D piecewise poly using adigatorGenInterp2pp and make it an input to your program and make your program dependent upon adigatorEvalInterp2pp. So, if you have a function
```
function y = myfunc(x,y,data)
    y = interp2(data.X,data.Y,data.Z,x,y,'spline')
end
```
you can now replace that with a function
```
function y = myfunc(x,data)
    y = interp2(data.pp2,x,y);
end
```
where data.pp2 can be generated apriori by the command data.pp2 = adigatorGenInterp2pp(data.X,data.Y,data.Z,'spline'). This should make your original file evaluate much faster.
* Analogous to the 1D ppval case, the adigatorEvalInterp2pp command may be used within a loop in the same way that ppval can be used (as described above). You shouldn't get any redundant calculations printed at the 2nd derivative level though.

V0.3.3 5/22/14 mjw
* User conditional statements containing short circuit and (&&) and/or shirt circuit or (||) will be replaced the non-short circuit versions (&) and (|), respectively, in the differentiated file - if the statements will not work with the non-short circuit versions then the file cannot be differentiated.
* Fixed an issue where if there is an if/elseif/else statement such that only one branch is always true, then there will not be a conditional statement printed to the derivative file, rather only the branch which is true.
* Overloaded the function atan2 and the corresponding code to cadabinaryarraymath.
* Fixed an issue when using the vectorized mode and assigning a vectorized object to an empty array.
* Fixed an issue with multiple cell/structure references on the same line within a loop.
Edited 5/26/14:
* A couple of the inverse trig functions (atan,acot,acsc,asec) had the wrong derivative rule coded - edited file cadaunaryarraymath.m to contain the proper rules.
Edited 5/27/14:
* There was an error being produced when the user had a single command spanning multiple lines (i.e. using "..." or building a matrix on multiple lines) and then commented them out, this should be fixed and the release was updated.

V0.3.2 5/15/14 mjw
* There was an issue using cell/structure references and assignments within loop statements - this should be fixed.
* Additionally, had a complaint about the README file, user manual, etc. being on the MATLAB path so I placed the utility functions (adigator, adigatorCreateDerivInput, etc.) within a /adigator/util/ folder and changed the startupadigator file to add this folder to the MATLAB path rather than /adigator/
V0.3.1
* Had everything coded for the '/' file separator - changed this so everything should work on windows machines
V0.3.0
* Changed the way the options work - options are now given directly to the adigator command
* Also made it so that adigator does not forcefully overwrite derivative files unless the overwrite option is set to true

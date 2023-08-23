Get Started with z-Tran
================
z-Tran is a MATLAB toolbox that Fei Tan and Jieran Wu developed for implementing our analytic policy function iteration (APFI) algorithm of the companion [paper](https://www.sciencedirect.com/science/article/abs/pii/S002205312100212X?via%3Dihub)&mdash;an algorithm that builds the algorithm for solving dynamic incomplete information models. All of the routines are coded in MATLAB as object-oriented programs. However, z-Tran does not require the user to have any prior knowledge of object-oriented programming. Indeed, it is quite intuitive to use as you will quickly realize when running our demo code. For any issue, suggestion or bug report, please send an email to [tanf [at] slu.edu](tanf@slu.edu).

Toolbox Structure
-----------------------------------
z-Tran consists of two classes defined under the "classes" folder: @varma (for **varma** value class) and @ztran (for **ztran** handle class).

1. The class definition file **varma.m** encapsulates the VARMA processes and the operations performed on them via a VARMA-type object.

2. The class definition file **ztran.m** encapsulates the model inputs (including the equilibrium equations, variable-innovation structure, and information structure) and its solving routine via a model-type object.

As is common in most numerical procedures, the algorithm convergence depends on many factors, e.g., the structure and stability of the model, the initial conjecture, the updating steps, the convergence criteria, etc. In the toolbox, we list a number of options that users can specify in the main interface **solve.m** to improve the numerical efficiency and stability. Alternatively, users can also construct their own code based on the APFI framework. Instead of using the canonical form and z-Tran interface, users call the individual routines under the "@varma" folder whenever is needed. This feature may be appealing in some complicated models when users need to monitor the iteration process to identify potential problems. See Online Appendix S3 for a practical user guide.

Replication Files
-----
To initialize z-Tran, set MATLAB directory to the toolbox's master folder and call **startup** at the command line, which will display the general information about z-Tran and add all relevant folders to MATLAB search path. The "examples" folder contains a growing collection of user files for solving the models considered in the literature:

* **Ex1_asset** - solves the asset pricing model in Section 2.

* **Ex2_dsge** - solves the DSGE model in Section 5.1.

* **Ex3_hank** - solves the HANK model in Section 5.2.

* **Ex4_rbc** - solves the RBC model in Section 5.3.

* **Ex5_confounding** - solves the beauty contest model in Section 2 of [Huo and Pedroni (2020)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3687529).

Acknowledgments
-----
We thank Zhao Han and Zheliang Zhu for testing z-Tran on various models and providing helpful user feedback.

Copyright Notice
-----
To avoid misinterpretation about z-Tran's contributors, we have updated the software's copyright notice (see [startup.m](https://github.com/econdojo/ztran/blob/main/startup.m) and [LICENSE](https://github.com/econdojo/ztran/blob/main/LICENSE)) under the BSD 3-Clause "New" or "Revised" License.

Suggested Citation
-----
Tan, Fei and Wu, Jieran, z-Tran: MATLAB Software for Solving Dynamic Incomplete Information Models (August 22, 2023). Available at GitHub: https://github.com/econdojo/ztran

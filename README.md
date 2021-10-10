Get Started with z-Tran
================
z-Tran is a MATLAB toolbox that implements the analytic policy function iteration (APFI) algorithm of [Han, Tan and Wu (2021)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3512320) for solving dynamic incomplete information models. All of the routines are coded in MATLAB as object-oriented programs. However, z-Tran does not require the user to have any prior knowledge of object-oriented programming. Indeed, it is quite intuitive to use as you will quickly realize when running our demo code. For any issue, suggestion or bug report, please send an email to [tanf [at] slu.edu](tanf@slu.edu)

Toolbox Structure
-----------------------------------
z-Tran consists of two classes defined under the "classes" folder: @varma (for **varma** class) and @ztran (for **ztran** class).

1. The class definition file **varma.m** encapsulates the VARMA processes and the operations performed on them via a VARMA-type object.

2. The class definition file **ztran.m** encapsulates the model inputs (including the equilibrium equations, variable-innovation structure, and information structure) and its solving routine via a model-type object.

As is common in most numerical procedures, the algorithm convergence depends on many factors, e.g., the structure and stability of the model, the initial conjecture, the updating steps, the convergence criteria, etc. In the toolbox, we list a number of options that users can specify in the main interface **solve.m** to improve the numerical efficiency and stability. Alternatively, users can also construct their own code based on the APFI framework. Instead of using the canonical form and z-Tran interface, users call the individual routines under the "@varma" folder whenever is needed. This feature may be appealing in some complicated models when users need to monitor the iteration process to identify potential problems. See Appendix C of [Han, Tan and Wu (2021)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3512320) for a practical user guide.

Replication Files
-----
To initialize z-Tran, set MATLAB directory to the toolbox's master folder and call **startup** at the command line, which will display the general information about z-Tran and add all relevant folders to MATLAB search path. The "examples" folder contains a growing collection of user files for solving the models considered in the literature:

* **Ex1_asset** - solves the asset pricing model in Section 2 of [Han, Tan and Wu (2021)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3512320).

* **Ex2_dsge** - solves the DSGE model in Section 5.1 of [Han, Tan and Wu (2021)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3512320).

* **Ex3_hank** - solves the HANK model in Section 5.2 of [Han, Tan and Wu (2021)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3512320).

* **Ex4_rbc** - solves the RBC model in Section 5.3 of [Han, Tan and Wu (2021)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3512320).

* **Ex5_confounding** - solves the beauty contest model in Section 2 of [Huo and Pedroni (2020)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3687529).

References
---------------------
Han, Z., F. Tan, and J. Wu (2021): “Analytic Policy Function Iteration,” Working Paper.

Huo, Z., and M. Pedroni (2020): “Dynamic Information Aggregation: Learning from the Past,” Working Paper.

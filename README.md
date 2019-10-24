# KBAL
Kernel Balancing project  

Package for implementation of kernel balancing. 

Investigators often use matching and weighting techniques to adjust for differences between treated and control groups on observed characteristics. These methods, however, require the user to choose what functions of the covariates must be balanced, and do not in general ensure equal multivariate densities of the treated and control groups. Treatment effect estimates made after adjustment by these methods are thus sensitive to specification choices, and are biased if any function of the covariates influencing the outcome has a different mean for the treated and control groups. This paper introduces kernel balancing, a method designed to reduce this bias without relying on specification searches or balance tests. The weights derived by kernel balancing (1) achieve approximate mean balance on a large class of smooth functions of the covariates, and (2) approximately equalize the multivariate densities of the treated and controls, when estimated a certain way. In two empirical applications, kernel balancing (1) accurately recovers the experimentally estimated effect of a job training program, and (2) finding that after controlling for observed differences, democracies are less likely to win counterinsurgencies, consistent with theoretical expectation but in contrast to previous findings.

See www.chadhazlett.com for details and paper. 

test
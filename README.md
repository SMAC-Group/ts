
An Introduction to Time Series Analysis with `R`
=========================================================

<center>
Stéphane Guerrier, Roberto Molinari & Haotian Xu
</center>
<br>


Introduction
------------

This text is **currently under development** to give readers a general overview of the main aspects that characterise time series and the most common tools used for their analysis. In particular, it is being used as a support for the STAT-463 class at Penn State University.

The book gives an overview of the practical settings in which time series analysis is of primary importance. Through the use of real-life examples and simulations, the book then proceeds to introducing the more formal definitions that are necessary to understand the problems behind time series analysis and the solutions that are available to those who are interested in performing such an analysis.

Starting from the definition of the basic time series models, the book provides the theoretical bases to understand the most important properties of time series among which the concept of stationarity and the different forms through which it is possible to characterise the main features of these processes.

Once the basic properties and characteristics are stated, the book proceeds to delivering an overview of the methods used to estimate the quantities of interest. From the empirical autocovariance functions to the MLE for time series model estimation, the reader will be lead through the technical details underlying the different estimators which allow them to better understand and interpret the outputs of an estimation procedure. Using the R software, the book will therefore allow the reader to carry out the main analyses to correctly address the complex features of time series.

Finally, a quick overview of more advanced topics in time series analysis is provided thereby covering heteroskedastic time series (GARCH), state-space models and, possibly, multivariate time series. The latter topics are also extremely common in the analysis of real-life time series analysis.


`R` and `RStudio`
-----------------

The statistical computing language `R` has become commonplace for many applications in industry, government, and academia. Having started as an open-source language to make available different statistics and analytical tools to researchers and the public, it steadily developed into one of the major software languages which not only allows to develop up-to-date, sound, and flexible analytical tools, but also to include these tools within a framework which is well-integrated with other important programming languages, communication, and version-control features. The latter is also possible thanks to the development of the `RStudio` interface which provides a pleasant and functional user-interface for `R` as well as an efficient Integrated Development Environment (IDE) in which different programming languages, web-applications and other important tools are available to the user. In order to illustrate the relationship between R & RStudio in statistical programming, one might think of a car analogy in which R would be like the engine and RStudio might be like leather seats. R is doing the work (for the most part), and RStudio generally is making you more comfortable while you use R.

Main References
---------------

This is not the first (or the last) book that has been written on time series analysis. Indeed, this can be seen as a book that brings together and reorganizes information and material from other sources structuring and tailoring it to a course in basic time series analysis. The main and excellent references (which are far from being an exhaustive review of literature) that can be used to have a more in-depth view of different aspects treated in this book are:

- Cochrane ([2005](http://econ.lse.ac.uk/staff/wdenhaan/teach/cochrane.pdf));
- Hamilton ([1994](http://doc1.lbfl.li/aca/FLMF037168.pdf));
- Shumway & Stoffer ([2010](http://db.ucsd.edu/static/TimeSeries.pdf)).

License
-------

You can redistribute it and/or modify this book under the terms of the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (CC BY-NC-SA) 4.0 License.

<a href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img src="/images/licence.png" align="left" width="200"/></a> <br><br><br>

References
==========

Cochrane, John H. 2005. *Time Series for Macroeconomics and Finance*.

Hamilton, James D. 1994. *Time series analysis*, 2nd ed. Princeton, NJ: Princeton University Press.

Shumway,  Robert H. & Stoffer, David S. 2010. *Time Series Analysis and Its Applications: With R Examples*. 3rd ed. New York, NY: Springer.

%% FDFIT TUTORIAL
% This tutorial will guide you through the basic functionality of the
% FDFIT package. We will analyze a sample dataset of double-stranded DNA
% force-extension curves, using the approach described in Ref. 1. Along the
% way, you will learn to load force-extension data into MATLAB, process and 
% visualize it, and fit the eWLC model to the data.

%% Before You Begin
% Make sure MATLAB's "Current Folder" is set to the folder containing the
% FDFIT package. That should be the folder containing a script called
% "init.m", as well as a bunch of subfolders.
% 
% Before you start, you should run the |init| script. This adds all the
% subfolders to MATLAB's 
% <matlab:web(fullfile(docroot,'matlab/search-path.html')) search path>, so
% MATLAB can find all the functions in the package.
% (Of course, it's also possible to do this yourself using MATLAB's user
% interface.)

init;


%% Loading the Sample Data
% The FDFIT package contains functions for loading force-extension data
% from various sources, including text files (CSV, ASCII, etc.). Let's load
% some sample data for double-stranded DNA, which has been included with
% this code.
%
% Note that all functions in this package are documented. See, for example,
% the documentation of the <matlab:doc('readasciifile') readasciifile>
% function used below.

fd = readasciifile('Doc/SampleData/Sample01.csv', 'HeaderLines', 6);


%% Plotting Data
% The |fd| variable we've just created contains our F,d (force-distance) 
% data. Let's have a look at it:

plotfd(fd);

%%
% <matlab:doc('plotfd') plotfd> supports various styles of plotting. For
% example, we can also plot the data as a pair of F,t and d,t graphs:

plotfd(fd, 'style', 'fdt');


%% Metadata
% The |fd| variable containing our data is a MATLAB _object_. In case
% you're not familiar with object-oriented programming: an object is
% similar to a |struct|, with a touch of added magic: besides just fields,
% the object can also contain associated functions.
%
% Let's have a closer look at the contents of |fd|:

disp(fd)

%%
% There isn't a lot to see, but your eye might be caught by the line
% "Metadata". As a matter of fact, |FdData| objects can keep track of
% simple metadata. For data from ASCII files, this isn't a lot, but we'll
% encounter some other examples where this is more useful.
%
% To show the metadata in its full glory, use:

disp(fd, 'full')

%%
% As you can see, the object contains, in its metadata, the original
% location the file was read from; as well as the header lines. This
% metadata was placed there by the |readasciifile| function.


%% Simple Data Corrections
% Now, while plotting the data, we have noticed that the distance values
% are off. Lambda DNA should have a contour length of 16.5 um, but judging
% from the F,d plot, it's more like 19.5 um. Looking back at the metadata,
% we can see where this is coming from: the distance values include the
% bead diameter. We'll have to subtract that.
%
% |FdData| objects contain some useful functions for these kinds of
% corrections. For example, to fix the bead diameter issue, we can use the
% <matlab:doc('FdData/shift') shift> function:

fd = fd.shift('d', -3.05);
figure;
plotfd(fd);

%%
% That looks more like it. Note that, if we now look at the |fd| variable
% again...

disp(fd)

%%
% The new |FdData| object that was returned by the |shift| function
% has remembered it's the shifted version of the original. Under
% "History", you see a line |"shift(d, -3.05)"| that wasn't there before.
%
% This is a useful little feature that should help with taking care of data
% provenance. All |FdData| functions that manipulate the data will add such
% a line to the history list, so you can retrace your steps later on. Some
% of these functions are:
%
% * <matlab:doc('FdData/shift') shift>: shift the data along a particular
%   axis, i.e., add a constant offset to a particular column of the data;
% * <matlab:doc('FdData/subset') subset>: take a subsest of the data, for
%   example, only the data between 0 and 30 pN, or only the data between 10
%   and 15 um.
% * <matlab:doc('FdData/scale') scale>: scale the data, i.e., multiply a
%   particular column by a scaling factor;
% * <matlab:doc('FdData/fragment') fragment>: takes a subrange of the data
%   ? like subindexing an array. Comparable to "subset", but uses data
%   indices instead of column values.

%%
% By the way, it's also possible to use a visual tool to select subsets of
% the data. Call the <matlab:doc('trimfd') trimfd> function to open a
% dialog that allows you to select a range of data using two draggable
% cursor lines:
%
%    fd_fragment = trimfd(fd);
%
% <<trimfd.png>>
%
% _The |trimfd| dialog window. Drag the red cursor lines in the right half
% of the screen to select a subset of the data._


%% Basic Fitting
% Of course, the main goal of the package is to fit F,d curves. Let's start
% by doing a simple Odijk eWLC fit to our sample curve.

f = fitfd(fd)

%%
% That's all: if you don't give any other options, this will fit the Odijk
% extensible worm-like chain [2,3], with a force offset, to all data up to 30
% pN. Let's plot the resulting fit:

plotfdfit(fd, f);

%%
% The <matlab:doc('fitfd') fitfd> function has quite a number of options
% (see the documentation for details). As one example, we could try fitting
% a slightly different fit model, and explicitly give the starting values
% for the fit parameters. Here, we'll use a model that includes an
% additional 'distance offset' parameter, so we allow for possible
% variations in bead diameter (i.e., horizontal shifts of the curve). In
% this case, we'll also have to provide the expected contour length (|Lc|),
% with respect to which the distance offset is measured.

f = fitfd(fd, 'model', 'odijk-d0-f0', 'startParams', [50 1500 0 0], 'Lc', 16.5)



%% Collections of Data
% Now you know some of the basic commands, we can turn to one of the more
% powerful features of the package: managing collections of data.
%
% We can load groups of |FdData| objects into an 'array', called an
% |FdDataCollection| (another type of object). Using these collections, we
% can easily explore the data, apply corrections on groups as a whole, and
% do global fits.
%
% Let's start by loading _all_ sample data curves into one of these
% collections. Oh, and while we're at it, let's not forget to correct for
% those bead diameters:

collection = readasciifolder('Doc/SampleData/*.csv', 'HeaderLines', 6, ...
    'beadDiameter', 3.05);

%%
% As before, we can use the |disp| function to peek inside:

disp(collection)

%%
% And we can also plot the data:

figure;
plotfd(collection);


%%
% Once we have data in a collection, there is also a simple user interface
% we can bring up to browse the data:
%
%    explorefd(collection);
% 
% <<FdExplorer.png>>
%
% In this window, there's a few nifty features:
%
% * You can plot multiple graphs simultaneously, by holding down the Ctrl key
%   (or the Cmd key on the Mac), and clicking in the data list on the
%   left-hand side.
% * The buttons in the top right corner allow you to visually trim the
%   selected data files (using the window with red cursors you've seen
%   before), and to remove files you don't like from the collection.
% * You can _tag_ data curves. Each curve in the collection can have one or
%   more associated tags, similar to the tags you may have encountered on
%   websites like Gmail. You can, for example, give all curves that were
%   measured in the presence of protein a tag |protein|. Later on, you can
%   then extract those protein curve using the command:
%
%     collection.getByTag('protein')
%
% This will give you a new collection containing just the F,d curves that
% have the |protein| tag.
%
% * The window also allows you to filter by tags, using the controls in the
%   bottom-left corner.

%%
% Collections are great, since they allow you to easily manipulate multiple
% datasets at once. For example, we can fit all the curvecs in one go, and
% then browse through the individual fit results:
%
%     fits = fitfd(collection);
%     explorefdfits(fits);
%
% In the _Cell Explorer_ that pops up, you can use the arrow, PgUp and PgDn
% keys on your keyboard to look at the fits we just made.

%%
% Most importantly, having data in a collection allows you to use the
% _global fitting_ technique described in the paper [1]. Simply run:

gfit = fitfdglobal(collection, 'Lc', 16.5)

%%
% If you want to inspect the fit results, use:
%
%     plotglobalfdfits(collection, gfit);



%% Twistable Worm-Like Chain Analysis
% The twistable worm-like chain (tWLC) analysis has been implemented as an
% object <matlab:doc('TwlcAnalysis') TwlcAnalysis>. To run a tWLC analysis
% on our sample data, follow these steps:

twa = TwlcAnalysis('data', collection, 'nBootstrapIter', 4)

%%
% This creates a TwlcAnalysis object, that will run all the steps described
% in the paper [1]. We'll decrease the number of bootstrap iterations for
% now, but for an actual analysis, you should probably keep this to the
% default setting of 100.
%
% For a precise description of all the parameters, I'll have to refer you
% to the paper [1], its Supplementary Information, as well as the
% <matlab:doc('TwlcAnalysis') TwlcAnalysis documentation>.
%
% As the output already indicates, we should now call the "analyze"
% function to run the actual analysis:

twa.analyze();

%%
% Plotting the fit:

twa.plotFit();

%%
% Usually, it's also a good idea to inspect the sweep results. You can do
% this using the <matlab:doc('TwlcAnalysis/plotSweeps') plotSweeps>
% function:

twa.plotSweeps();


%% REFERENCES
%
% # O.D. Broekmans, G.A. King, G.J. Stephens, G.J.L. Wuite, 
%   _DNA Twist Changes With Magnesium(2+) Concentration_,
%   arXiv:1412.8753 [physics.bio-ph].
% # T. Odijk, Stiff Chains and Filaments under Tension, Macromolecules
%   28, 7016-7018 (1995).
% # M. D. Wang, H. Yin, R. Landick, J. Gelles, S. M. Block, Stretching
%   DNA with optical tweezers., Biophysical journal 72, 1335-46 (1997).

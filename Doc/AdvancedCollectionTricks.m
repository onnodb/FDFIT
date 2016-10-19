%% ADVANCED COLLECTION TRICKS
% For more experienced MATLAB users, the |FdDataCollection| object harbors
% some nice time-saving tricks. In this small tutorial, I'll show you some
% of the most useful ones.
%
% During this tutorial, I'll assume you know what a
% <matlab:web(fullfile(docroot,'matlab/matlab_prog/overview-1.html'))
% function handle> is. You'll also need to have some basic familiarity with
% <matlab:web(fullfile(docroot,'matlab/matlab_prog/anonymous-functions.html'))
% anonymous functions>.

%% Loading the sample data
% Let's start by creating an |FdDataCollection| with the sample data:

myCollection = readasciifolder('Doc/SampleData/*.csv', 'HeaderLines', 6);

%%

disp(myCollection);

%%
% A plot of the sample data:

plotfd(myCollection);


%% Using Tags
% As mentioned in the tutorial, |FdData| objects can be associated with
% _tags_ (or 'categories'; you've probably encountered them on the internet
% before, e.g., Gmail). For the sake of this example, we will add some tags
% ourselves:

myCollection.items{2}.addTag('looksnice');
myCollection.items{4}.addTag('looksnice');
myCollection.items{5}.addTag('looksnice', 'favorite');
myCollection.items{6}.addTag('favorite');
myCollection.items{10}.addTag('shortOSPlateau');
disp(myCollection);

%%
% The table produced by the <matlab:doc('FdDataCollection/disp') disp>
% function contains an overview of the tags we just added.
%
% We can now easily retrieve curves from the collection that have a
% particular tag:

nice = myCollection.getByTag('looksnice');

%%

disp(nice);

%%
% What's important to realize here, is that the objects (F,d curves) in our
% new |nice| collection, are still the same ones as in the original
% |myCollection|. They're just part of multiple |FdDataCollections| _at the
% same time_. (However, a collection cannot contain the same item more than
% once). We can see the effect of this if we add a tag to one of the items
% in |nice|:

disp(nice);

%%

disp(myCollection)

%%

nice.items{2}.addTag('supercalifragilistic');
disp(nice);

%%

disp(myCollection);

%%
% As you can see, the F,d curve was updated in both collections, _because
% they're really the same object_. If you want to read more about this, see
% the MATLAB documentation on
% <matlab:web(fullfile(docroot,'matlab/ref/handle.html')) handle objects>.
%
% Leaving that little detour, let's summarize the main point: you can use
% tags to organize your data. Instead of keeping, for example, different
% magnesium concentrations in different collections, you can just stuff all
% the data into one bit |FdDataCollection|, and then easily retrieve data
% for a particular concentration using the
% <matlab:doc('FdDataCollection/getByTag') getByTag> function.


%% Combining Collections
% The |FdDataCollection| also contains a number of functions for combining
% collections in various ways. Say, for example, that we'd now like to find
% all F,d curves that have _either_ the |looksnice| _or_ the |favorite|
% tag. How do we do that?
%
% Simple: we get collections for each of those tags separately, and then
% just 'add them up'. In set theory terms: we take the
% <matlab:doc('FdDataCollection/union') union>.

nice = myCollection.getByTag('looksnice');
favs = myCollection.getByTag('favorite');
niceAndFavs = nice.union(favs);
disp(niceAndFavs);

%%
% Or, in a single line:
%
%    disp(myCollection.getByTag('looksnice').union(myCollection.getByTag('favorite')));
%
% You'll find more functions like this, including
% <matlab:doc('FdDataCollection/intersect') intersect>,
% <matlab:doc('FdDataCollection/subtract') subtract>,
% <matlab:doc('FdDataCollection/isequal') isequal>, and
% <matlab:doc('FdDataCollection/isempty') isempty>.


%% Applying Operations to Collections
% Finally, the last time-saving trick involves the ability to easily apply
% operations to all curves in a collection.
%
% Say we'd like to get subsets of the F,d curves: only the data below 30
% pN. We could get this by writing a |for| loop:

dataBelow30pN = FdDataCollection();
for i = 1:myCollection.length
    dataBelow30pN.add(myCollection.items{i}.subset('f', [-Inf 30]));
end

%%
% That works perfectly fine, but you might enjoy using this syntax more:

dataBelow30pN = myCollection.map(@(fd) fd.subset('f', [-Inf 30]));

%%
% Or, if you'd like to replace all items in |myCollection|, instead of
% getting a new collection with subsetted curves:

myCollection.applyToAll(@(fd) fd.subset('f', [-Inf 30]));

%%
% To complement these two functions, there are also
% <matlab:doc('FdDataCollection/filter') filter>, to get a collection
% containing only items for which a particular function returns |true|; and
% <matlab:doc('FdDataCollection/forAll') for>, for running a function that
% doesn't return anything on all items in the collection.


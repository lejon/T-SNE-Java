T-SNE-Java
==========


Pure Java implementation of Van Der Maaten and Hinton's t-sne clustering algorithm.

This is a simple port of the Python version of T-SNE.

There are two versions, one the fastest, uses the EJML library for fast matrix operations, this version is called FastTSne.
The second version is a more direct port, but slower with call-by-value semantics, this version is called SimpleTSne.
Both versions implement the TSne interface.

There is still room for improvements, cleanup is the first that comes to mind. :)

The easiest way to play around with it right now, still involves installin Apache Maven. Hopefully will get to doing a
proper release at GitHub soon, but for now.

Install Apache Maven, then:

	mvn package
	
To run a demo:

	java -cp target/tsne-X.X.X-SNAPSHOT.jar com.jujutsu.tsne.TSneDemo


This is still a "development" version, i.e it is currently not too useful as a stand alone tool.

To see graps generated with this implementation, [Klick here](http://lejon.github.io/TSneJava/)

Enjoy!
-Leif
  
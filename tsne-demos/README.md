T-SNE-Java Demos
================


Pure Java implementation of Van Der Maaten and Hinton's t-sne clustering algorithm.

The core packages are in "tsne-core" from the toplevel dir. 
	
To run a demo:

	java -cp target/tsne-demos-0.0.1-SNAPSHOT.jar com.jujutsu.tsne.demos.TSneDemo mnist2500_X.txt mnist2500_labels.txt
	
The two example files are in 'src/main/resources/datasets'.

This is still a "development" version, i.e it is currently not too useful as a stand alone tool.

To see graphs generated with this implementation, [Klick here](http://lejon.github.io/TSneJava/)

Enjoy!
-Leif
  

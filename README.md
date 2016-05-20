[![Build Status](https://travis-ci.org/lejon/T-SNE-Java.svg?branch=master)](https://travis-ci.org/lejon/T-SNE-Java)

T-SNE-Java
==========


Pure Java implementation of Van Der Maaten and Hinton's t-sne clustering algorithm.

This project has been refacored into two separate Maven projects, one for the core tsne and one for the demos.

T-SNE-Java is now in version 1.0 (tsne-core is 1.0.1)  It has been in development version long enough.

Basic command line usage
------------------------

If you just want to use TSne as a command line tool, you should use TSneCsv which is by far the most flexible alternative.

You must then first build and install 'tsne' and  'tsne-demos' (mvn install).
Then use the tsne-demos JAR just build. 

You can also download the pre-build binary JAR from the release page.

Examples:

Run TSne on file without headers and no labels.
```shell
java -jar target/tsne-demos-1.0.jar -nohdr -nolbls src/main/resources/datasets/iris_X.txt 
```
Run TSne on CSV file with headers and label column nr. 5.
```shell
java -jar target/tsne-demos-1.0.jar --lblcolno 5 src/main/resources/datasets/iris.csv
```
Run TSne on file without headers and no labels but supply a separate label file (with the same ordering as the data file).
```shell
java -jar target/tsne-demos-1.0.jar --nohdr --nolbls --label_file=src/main/resources/datasets/iris_X_labels.txt src/main/resources/datasets/iris_X.txt
```

To see graphs generated with this implementation, [Klick here](http://lejon.github.io/TSneJava/)

For some tips working with t-sne [Klick here] (http://lejon.github.io)

Basic code usage: 
-----------------

```java
import com.jujutsu.tsne.FastTSne;
import com.jujutsu.tsne.MatrixOps;
import com.jujutsu.tsne.TSne;

public class TSneTest {
  public static void main(String [] args) {
    int initial_dims = 55;
    double perplexity = 20.0;
    double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/resources/datasets/mnist2500_X.txt"), ",");
    System.out.println(MatrixOps.doubleArrayToPrintString(X, ", ", 50,10));
    TSne tsne = new FastTSne();
    double [][] Y = tsne.tsne(X, 2, initial_dims, perplexity);   
    
    // Plot Y or save Y to file and plot with some other tool such as for instance R
    
  }
}

```

Version
-------
Demo: 1.0
Core: 1.0.1

Enjoy!
-Leif
  

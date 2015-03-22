T-SNE-Java
==========


Pure Java implementation of Van Der Maaten and Hinton's t-sne clustering algorithm.

This project has been refacored into two separate Maven projects, one for the core tsne and one for the demos.

This is still a "development" version, i.e it is currently not too useful as a stand alone tool.

Basic usage: 
  
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

To see graps generated with this implementation, [Klick here](http://lejon.github.io/TSneJava/)

For some tips working with t-sne [Klick here] (http://lejon.github.io)

Enjoy!
-Leif
  

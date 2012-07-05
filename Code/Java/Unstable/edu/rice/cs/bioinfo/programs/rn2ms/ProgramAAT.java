package edu.rice.cs.bioinfo.programs.rn2ms;
import edu.rice.cs.bioinfo.library.programming.Proc1;
import org.junit.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/29/12
 * Time: 2:40 PM
 * To change this template use File | Settings | File Templates.
 */
public class ProgramAAT {

   private StringBuffer outBuffer;

   private StringBuffer errorBuffer;

   private Proc1<String> out;

   private Proc1<String> error;

    @Before public void setUp()
    {
        outBuffer   = new StringBuffer();
        errorBuffer = new StringBuffer();

        out = new Proc1<String>() {
            public void execute(String input) {
                outBuffer.append(input);
            }
        };

        error = new Proc1<String>() {
            public void execute(String input) {
                errorBuffer.append(input);
            }
        };


    }

    @Test
    public void testRun1()
    {
        String[] args = new String[] { "10", "((u#1:0::.3,A:2)x:1,(((B:1,C:1)v:1)u#1:0::.7,D:2)y:1)R;" };
        Program.run(args, out, error);
        String result = outBuffer.toString();
        Assert.assertEquals("Node to population:\n" +
                "{D=4, A=1, B=2, C=3}\n" +
                "MS command:\n" +
                "ms 4 10 -T -I 4 1 1 1 1 -ej 1 3 2 -ej 2 1 2 -es 2 2 0.3 -ej 2 4 5 -ej 3 2 5 ", result);

    }
}

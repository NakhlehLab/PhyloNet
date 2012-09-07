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
        String[] args = new String[] { "10", "((u#1:0::.3,A:2)x:1,(((B:1,C:1)v:1)u#1:0::.7,D:2)y:1)R;", "0" };
        Program.run(args, out, error);
        String result = outBuffer.toString();
        Assert.assertEquals("Node to population:\n" +
                "{D=4, A=1, B=2, C=3}\n" +
                "MS command:\n" +
                "ms 4 10 -T -I 4 1 1 1 1 -ej 1 3 2 -ej 2 1 2 -es 2 2 0.3 -ej 2 4 5 -ej 3 2 5 ", result);

    }

    @Test
    public void testRun2()
    {
        String[] args = new String[] { "10", "(_6#H1:0.15886101::0.08886512207602625945668251006281934678554534912109375,(((5:0.26578128,(4:0.02449977,3:0.02449977)_4:0.24128151)_3:0.15886101)_6#H1:0::0.91113487792397374054331748993718065321445465087890625,(2:0.12419069,1:0.12419069)_2:0.30045161)_5:0.15886101)_1;", ".1" };
        Program.run(args, out, error);
        String result = outBuffer.toString();
    }
}

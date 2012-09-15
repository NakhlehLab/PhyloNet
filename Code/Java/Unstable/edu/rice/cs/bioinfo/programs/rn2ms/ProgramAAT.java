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
                "ms 4 10 -T -I 4 1 1 1 1 -ej 0.5 3 2 -ej 1 1 2 -es 1 2 0.3 -ej 1 4 5 -ej 1.5 2 5 ", result);

    }

    @Test
    public void testRun2()
    {
        String[] args = new String[] { "10", "((_9#H3:0.11681459625::0.9445572341233876567656579936738125979900360107421875,(_5#H1:0.049845455::0.22415176386098456529083478017128072679042816162109375,2:0.09969091)_2:0.22965631)_8:0.22965631,((_7#H2:0.18760989625::0.71198742555209182381048549359547905623912811279296875)_9#H3:0.18760989625::0.0554427658766123432343420063261874020099639892578125,((((1:0.0249227275)_7#H2:0.0249227275::0.28801257444790817618951450640452094376087188720703125)_5#H1:0.070795300::0.77584823613901543470916521982871927320957183837890625,3:0.120640755)_4:0.120640755,4:0.24128151)_3:0.15886101)_6:0.15886101)_1;", ".1" };
        Program.run(args, out, error);
        String result = outBuffer.toString();
    }
}

package edu.rice.cs.bioinfo.programs.networksearchgen;

import edu.rice.cs.bioinfo.library.programming.Proc1;
import edu.rice.cs.bioinfo.programs.rn2ms.*;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/7/12
 * Time: 2:19 PM
 * To change this template use File | Settings | File Templates.
 */
public class ProgramAAT
{
    private StringBuffer outBuffer;

    private Proc1<String> out;

    @Before
    public void setUp()
    {
        outBuffer   = new StringBuffer();

        out = new Proc1<String>() {
            public void execute(String input) {
                outBuffer.append(input);
            }
        };
    }

    @Test
    public void testRun1()
    {
        String[] args = new String[] { "ms", "10", "((4:2.41281510,3:2.41281510)_3:3.17722020,(2:0.99690910,1:0.99690910)_2:4.59312620)_1;", ".0001" };
        edu.rice.cs.bioinfo.programs.networksearchgen.Program.run(args, out);
        String result = outBuffer.toString();
        Assert.assertEquals("{3=2, 2=3, 1=4, 4=1}\n" +
                "ms 4 10 -T -I 4 1 1 1 1 -ej 0.4984545499999999407236828119494020938873291015625 4 3 -ej 1.206407550000000217238493860349990427494049072265625 2 1 -ej 2.7950176500000001311718733632005751132965087890625 1 3".trim(),
                result.trim());

    }

    @Test
    public void testRun2()
    {
        String[] args = new String[] { "st", "((3:1.066503,4:1.066503):5.188327,(1:2.534810,2:2.534810):3.720020);//((3:1.066503,4:1.066503):5.188327,(1:2.534810,2:2.534810):3.720020);//((3:1.004552,4:1.004552):4.696378,(1:2.496120,2:2.496120):3.204810);//((3:1.806840,4:1.806840):4.616031,(1:2.906498,2:2.906498):3.516373);//((3:1.171307,4:1.171307):6.171578,(1:2.769297,2:2.769297):4.573588);"};
        edu.rice.cs.bioinfo.programs.networksearchgen.Program.run(args, out);
        String result = outBuffer.toString();
        Assert.assertEquals("((2:5.069620,1:5.069620)_2:7.440040,(4:2.133006,3:2.133006)_1:10.376654)_0;\n" +
                "((2:5.069620,1:5.069620)_5:7.440040,(4:2.133006,3:2.133006)_4:10.376654)_3;\n" +
                "((2:4.992240,1:4.992240)_8:6.409620,(4:2.009104,3:2.009104)_7:9.392756)_6;\n" +
                "((2:5.812996,1:5.812996)_11:7.032746,(4:3.613680,3:3.613680)_10:9.232062)_9;\n" +
                "((2:5.538594,1:5.538594)_14:9.147176,(4:2.342614,3:2.342614)_13:12.343156)_12;", result.trim());

    }





}

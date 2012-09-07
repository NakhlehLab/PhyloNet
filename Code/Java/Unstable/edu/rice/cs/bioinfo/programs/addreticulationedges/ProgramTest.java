package edu.rice.cs.bioinfo.programs.addreticulationedges;
import org.junit.*;

import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/13/12
 * Time: 4:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class ProgramTest
{
    @Test
    public void testAddEdges()
    {
        String result = new Program().addEdges("(A:10,B:20)R;", 1, 1.0, new Random(13));
        Assert.assertEquals("(_2#H1:5::0.32885982758721066243623454283806495368480682373046875,((B:15)_2#H1:0::0.67114017241278933756376545716193504631519317626953125,A:5)_1:5)R;", result);

        result = new Program().addEdges("((C:5,D:5)A:10,(E:2,F:2)B:20)R;", 1, 1.0, new Random(15));
        Assert.assertEquals("(_2#H1:12.5::0.7351050399591991801884205415262840688228607177734375,((((F:2,E:2)B:7.5)_2#H1:0::0.2648949600408008198115794584737159311771392822265625,C:2.5)_1:2.5,D:5)A:10)R;", result);

        result = new Program().addEdges("((C:5,D:5)A:10,(E:2,F:2)B:20)R;", 1, 2.0, new Random(15));
        Assert.assertEquals("(_2#H1:25.0::0.7351050399591991801884205415262840688228607177734375,((((F:4,E:4)B:15.0)_2#H1:0::0.2648949600408008198115794584737159311771392822265625,C:5.0)_1:5.0,D:10)A:20)R;", result);

        result = new Program().addEdges("((C:5,D:5)A:10,(E:2,F:2)B:20)R;", 2, 1.0, new Random(15));
        Assert.assertEquals("(_2#H1:12.5::0.7329735365637322086485028194147162139415740966796875,(_4#H2:3.75::0.6034557804828082350212525852839462459087371826171875,((((D:1.25)_4#H2:0::0.3965442195171917649787474147160537540912628173828125,(F:2,E:2)B:6.25)_3:1.25)_2#H1:0::0.2670264634362677913514971805852837860584259033203125,C:2.5)_1:2.5)A:10)R;", result);

        result = new Program().addEdges("(1:0.45931262,(2:0.14159060,3:0.14159060):0.31772202);", 1, 1.0, new Random(13));
    }
}

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
        Assert.assertEquals("(_2#H1:5.0,((B:15.0)_2#H1:0.0,A:5.0)_1:5.0)R;", result);

        result = new Program().addEdges("((C:5,D:5)A:10,(E:2,F:2)B:20)R;", 1, 1.0, new Random(15));
        Assert.assertEquals("(_2#H1:12.5,((((F:2.0,E:2.0)B:7.5)_2#H1:0.0,C:2.5)_1:2.5,D:5.0)A:10.0)R;", result);

        result = new Program().addEdges("((C:5,D:5)A:10,(E:2,F:2)B:20)R;", 1, 2.0, new Random(15));
        Assert.assertEquals("(_2#H1:25.0,((((F:4.0,E:4.0)B:15.0)_2#H1:0.0,C:5.0)_1:5.0,D:10.0)A:20.0)R;", result);

        result = new Program().addEdges("((C:5,D:5)A:10,(E:2,F:2)B:20)R;", 2, 1.0, new Random(15));
        Assert.assertEquals("(_2#H1:12.5,(_4#H2:3.75,((((D:1.25)_4#H2:0.0,(F:2.0,E:2.0)B:6.25)_3:1.25)_2#H1:0.0,C:2.5)_1:2.5)A:10.0)R;", result);

        result = new Program().addEdges("(1:0.45931262,(2:0.14159060,3:0.14159060):0.31772202);", 1, 1.0, new Random(13));
    }
}

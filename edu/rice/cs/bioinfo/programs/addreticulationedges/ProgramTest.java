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
        Assert.assertEquals("((B:10)i2#H1:10::0.32885982758721066243623454283806495368480682373046875,(i2#H1:5::0.67114017241278933756376545716193504631519317626953125,A:5)i1:5)R;", result);

        result = new Program().addEdges("((C:5,D:5)A:10,(E:2,F:2)B:20)R;", 1, 1.0, new Random(15));
        Assert.assertEquals("(((E:1)i2#H1:16::0.38062203043353715070651333007845096290111541748046875,(D:5,C:5)A:5)i1:5,(i2#H1:1::0.61937796956646284929348666992154903709888458251953125,F:2)B:20)R;", result);

        result = new Program().addEdges("((C:5,D:5)A:10,(E:2,F:2)B:20)R;", 1, 2.0, new Random(15));
        Assert.assertEquals("(((E:2)i2#H1:32::0.38062203043353715070651333007845096290111541748046875,(D:10,C:10)A:10)i1:10,(i2#H1:2::0.61937796956646284929348666992154903709888458251953125,F:4)B:40)R;", result);

        result = new Program().addEdges("((C:5,D:5)A:10,(E:2,F:2)B:20)R;", 2, 1.0, new Random(15));
        Assert.assertEquals("((((F:1)i4#H2:8::0.1689738364905790657388706677011214196681976318359375,(E:1)i2#H1:8::0.19175955892193741192386369220912456512451171875)i3:8,(D:5,C:5)A:5)i1:5,(i4#H2:1::0.8310261635094209342611293322988785803318023681640625,i2#H1:1::0.80824044107806258807613630779087543487548828125)B:20)R;", result);

    }
}

package edu.rice.cs.bioinfo.programs.richnewick2hybridsimnewick;

import junit.framework.Assert;
import org.junit.Test;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/26/12
 * Time: 5:07 PM
 * To change this template use File | Settings | File Templates.
 */
public class ProgramTest
{
    @Test
    public void testConvertString1()
    {
        String result = Program.convertString("(((3:0.120640755)i5#H1:0.208706465::0.5489856822445762229989441038924269378185272216796875,(2:0.09969091,1:0.09969091)i2:0.22965631)i4:0.22965631,(i5#H1:0.120640755::0.4510143177554237770010558961075730621814727783203125,4:0.24128151)i3:0.31772202)i1;", "500");
        Assert.assertEquals("", result);
    }
}

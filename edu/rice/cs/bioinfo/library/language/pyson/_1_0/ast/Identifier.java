package  edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/31/11
 * Time: 1:44 PM
 * To change this template use File | Settings | File Templates.
 */
public class Identifier extends PySONNodeLineAndCol
{
    public final String Content;

    public Identifier(String content, int line, int col)
    {
        super(line, col);
        Content = content;
    }
}

package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/14/11
 * Time: 4:12 PM
 * To change this template use File | Settings | File Templates.
 */
public class Text implements AbstractSyntaxNode
{
    public final String Content;

    public final boolean OriginallyQuoted;

    public Text(String content, boolean originallyQuoted)
    {
        Content = content;
        OriginallyQuoted = originallyQuoted;
    }
}

package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

import com.sun.org.apache.bcel.internal.classfile.LineNumber;

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

    public final int LineNumberStart;

    public final int ColumnNumberStart;

    public Text(String content, int lineNumberStart, int columnNumberStart, boolean originallyQuoted)
    {
        Content = content;
        LineNumberStart = lineNumberStart;
        ColumnNumberStart = columnNumberStart;
        OriginallyQuoted = originallyQuoted;
    }
}

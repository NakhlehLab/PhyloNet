package edu.rice.bioinfo.library.language.pyson._1_0.ast;

import java.nio.ReadOnlyBufferException;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 5:30 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class PySONNodeLineAndCol implements PySONNode {

    public final int Line, Col;

      PySONNodeLineAndCol(int line, int col)
      {
         Line = line;
         Col = col;
      }
}

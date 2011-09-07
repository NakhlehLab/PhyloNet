/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/1/11
 * Time: 11:27 AM
 * To change this template use File | Settings | File Templates.
 */

package edu.rice.bioinfo.library.language.parsing;

public interface CoordinateParseError {

    public String getMessage();

    public int getLineNumber();

    public int getColumnNumber();


}

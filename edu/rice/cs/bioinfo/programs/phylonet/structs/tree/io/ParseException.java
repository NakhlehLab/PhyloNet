package edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io;

/**
 * This exception is thrown by readers when a tree is incorrectedly formatted.
 * 
 * @author Derek Ruths
 */
public class ParseException extends Exception {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public ParseException(String msg) {
		super(msg);
	}
}

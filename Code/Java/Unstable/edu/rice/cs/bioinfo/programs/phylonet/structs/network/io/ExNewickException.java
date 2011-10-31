package edu.rice.cs.bioinfo.programs.phylonet.structs.network.io;

/**
 * This class is used in parsing eNewick-format network input.
 */
public class ExNewickException extends RuntimeException {
	// Constructors
	/**
	 * Instantiates a new <code>ExNewickException</code> with the given
	 * message.
	 * 
	 * @param msg	the specific reason for this exception's existence
	 */
	public ExNewickException(String msg)
	{
		super(msg);
	}
	
	// Data members
	final static long serialVersionUID = 1;
}

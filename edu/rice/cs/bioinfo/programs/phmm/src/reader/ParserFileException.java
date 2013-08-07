package reader;

public class ParserFileException extends RuntimeException {

    /**
     *	Exception is thrown any time an incorrect read is done while parsing sequencing file.
     */
    private static final long serialVersionUID = 1L;

    //Parameterless Constructor
      public ParserFileException() {}

      //Constructor that accepts a message
      public ParserFileException(String message)
      {
         super(message);
      }


}

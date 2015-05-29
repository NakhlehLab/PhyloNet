// $ANTLR 3.4 D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g 2015-05-29 14:53:08

package edu.rice.cs.bioinfo.library.language.pyson._1_0.parsers.antlr.ast;
import java.util.LinkedList;


import org.antlr.runtime.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;

@SuppressWarnings({"all", "warnings", "unchecked"})
public class PySONLexer extends Lexer {
    public static final int EOF=-1;
    public static final int T__34=34;
    public static final int T__35=35;
    public static final int T__36=36;
    public static final int T__37=37;
    public static final int T__38=38;
    public static final int T__39=39;
    public static final int T__40=40;
    public static final int T__41=41;
    public static final int BEGIN=4;
    public static final int DATA=5;
    public static final int DATATYPE=6;
    public static final int DEFAULT_INDICATOR=7;
    public static final int DIMENSIONS=8;
    public static final int ELSE=9;
    public static final int END=10;
    public static final int FORMAT=11;
    public static final int GAP=12;
    public static final int ID=13;
    public static final int ID_SET=14;
    public static final int MATRIX=15;
    public static final int MISSING=16;
    public static final int MORPHDATA=17;
    public static final int NCHAR=18;
    public static final int NESTED_ML_COMMENT=19;
    public static final int NETWORK=20;
    public static final int NETWORKS=21;
    public static final int NTAX=22;
    public static final int PHYLONET=23;
    public static final int QUOTE=24;
    public static final int RN_LS_NONCOMMENT=25;
    public static final int START=26;
    public static final int SYMBOLS=27;
    public static final int TAXON_SET_LIST=28;
    public static final int TRANSLATE=29;
    public static final int TREE=30;
    public static final int TREES=31;
    public static final int UTREE=32;
    public static final int WS=33;

    public class ErrorWrapper
    {
    	public final String Message;
    	public final int Line,Col;
    	
    	ErrorWrapper(String message, int line, int col)
    	{
    		Message = message;
    		Line = line;
    		Col = col;
    	}
    }

    private List<ErrorWrapper> errors = new LinkedList<ErrorWrapper>();

    @Override   
    public void displayRecognitionError(String[] tokenNames, RecognitionException e) 
    {
            errors.add(new ErrorWrapper(getErrorMessage(e, tokenNames), e.line, e.c));
    }
        
    public List<ErrorWrapper> getErrors() 
    {
            return errors;
    }


    // delegates
    // delegators
    public Lexer[] getDelegates() {
        return new Lexer[] {};
    }

    public PySONLexer() {} 
    public PySONLexer(CharStream input) {
        this(input, new RecognizerSharedState());
    }
    public PySONLexer(CharStream input, RecognizerSharedState state) {
        super(input,state);
    }
    public String getGrammarFileName() { return "D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g"; }

    // $ANTLR start "T__34"
    public final void mT__34() throws RecognitionException {
        try {
            int _type = T__34;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:34:7: ( '(' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:34:9: '('
            {
            match('('); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__34"

    // $ANTLR start "T__35"
    public final void mT__35() throws RecognitionException {
        try {
            int _type = T__35;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:35:7: ( ')' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:35:9: ')'
            {
            match(')'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__35"

    // $ANTLR start "T__36"
    public final void mT__36() throws RecognitionException {
        try {
            int _type = T__36;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:36:7: ( ',' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:36:9: ','
            {
            match(','); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__36"

    // $ANTLR start "T__37"
    public final void mT__37() throws RecognitionException {
        try {
            int _type = T__37;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:37:7: ( ':' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:37:9: ':'
            {
            match(':'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__37"

    // $ANTLR start "T__38"
    public final void mT__38() throws RecognitionException {
        try {
            int _type = T__38;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:38:7: ( ';' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:38:9: ';'
            {
            match(';'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__38"

    // $ANTLR start "T__39"
    public final void mT__39() throws RecognitionException {
        try {
            int _type = T__39;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:39:7: ( '<' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:39:9: '<'
            {
            match('<'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__39"

    // $ANTLR start "T__40"
    public final void mT__40() throws RecognitionException {
        try {
            int _type = T__40;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:40:7: ( '=' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:40:9: '='
            {
            match('='); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__40"

    // $ANTLR start "T__41"
    public final void mT__41() throws RecognitionException {
        try {
            int _type = T__41;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:41:7: ( '>' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:41:9: '>'
            {
            match('>'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__41"

    // $ANTLR start "DEFAULT_INDICATOR"
    public final void mDEFAULT_INDICATOR() throws RecognitionException {
        try {
            int _type = DEFAULT_INDICATOR;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:170:2: ( '*' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:170:4: '*'
            {
            match('*'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "DEFAULT_INDICATOR"

    // $ANTLR start "PHYLONET"
    public final void mPHYLONET() throws RecognitionException {
        try {
            int _type = PHYLONET;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:172:9: ( ( 'P' | 'p' ) ( 'H' | 'h' ) ( 'Y' | 'y' ) ( 'L' | 'l' ) ( 'O' | 'o' ) ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:172:11: ( 'P' | 'p' ) ( 'H' | 'h' ) ( 'Y' | 'y' ) ( 'L' | 'l' ) ( 'O' | 'o' ) ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' )
            {
            if ( input.LA(1)=='P'||input.LA(1)=='p' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='H'||input.LA(1)=='h' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='Y'||input.LA(1)=='y' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='L'||input.LA(1)=='l' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='O'||input.LA(1)=='o' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "PHYLONET"

    // $ANTLR start "DATA"
    public final void mDATA() throws RecognitionException {
        try {
            int _type = DATA;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:176:6: ( ( 'D' | 'd' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'A' | 'a' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:176:8: ( 'D' | 'd' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'A' | 'a' )
            {
            if ( input.LA(1)=='D'||input.LA(1)=='d' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='A'||input.LA(1)=='a' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='A'||input.LA(1)=='a' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "DATA"

    // $ANTLR start "MORPHDATA"
    public final void mMORPHDATA() throws RecognitionException {
        try {
            int _type = MORPHDATA;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:179:2: ( ( 'M' | 'm' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'P' | 'p' ) ( 'H' | 'h' ) ( 'D' | 'd' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'A' | 'a' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:179:4: ( 'M' | 'm' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'P' | 'p' ) ( 'H' | 'h' ) ( 'D' | 'd' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'A' | 'a' )
            {
            if ( input.LA(1)=='M'||input.LA(1)=='m' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='O'||input.LA(1)=='o' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='R'||input.LA(1)=='r' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='P'||input.LA(1)=='p' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='H'||input.LA(1)=='h' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='D'||input.LA(1)=='d' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='A'||input.LA(1)=='a' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='A'||input.LA(1)=='a' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "MORPHDATA"

    // $ANTLR start "DIMENSIONS"
    public final void mDIMENSIONS() throws RecognitionException {
        try {
            int _type = DIMENSIONS;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:182:2: ( ( 'D' | 'd' ) ( 'I' | 'i' ) ( 'M' | 'm' ) ( 'E' | 'e' ) ( 'N' | 'n' ) ( 'S' | 's' ) ( 'I' | 'i' ) ( 'O' | 'o' ) ( 'N' | 'n' ) ( 'S' | 's' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:182:4: ( 'D' | 'd' ) ( 'I' | 'i' ) ( 'M' | 'm' ) ( 'E' | 'e' ) ( 'N' | 'n' ) ( 'S' | 's' ) ( 'I' | 'i' ) ( 'O' | 'o' ) ( 'N' | 'n' ) ( 'S' | 's' )
            {
            if ( input.LA(1)=='D'||input.LA(1)=='d' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='I'||input.LA(1)=='i' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='M'||input.LA(1)=='m' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='S'||input.LA(1)=='s' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='I'||input.LA(1)=='i' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='O'||input.LA(1)=='o' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='S'||input.LA(1)=='s' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "DIMENSIONS"

    // $ANTLR start "TREE"
    public final void mTREE() throws RecognitionException {
        try {
            int _type = TREE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:185:7: ( ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:185:10: ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' )
            {
            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='R'||input.LA(1)=='r' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "TREE"

    // $ANTLR start "UTREE"
    public final void mUTREE() throws RecognitionException {
        try {
            int _type = UTREE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:187:8: ( ( 'U' | 'u' ) ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:187:11: ( 'U' | 'u' ) ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' )
            {
            if ( input.LA(1)=='U'||input.LA(1)=='u' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='R'||input.LA(1)=='r' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "UTREE"

    // $ANTLR start "NETWORK"
    public final void mNETWORK() throws RecognitionException {
        try {
            int _type = NETWORK;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:189:9: ( ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:189:11: ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' )
            {
            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='W'||input.LA(1)=='w' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='O'||input.LA(1)=='o' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='R'||input.LA(1)=='r' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='K'||input.LA(1)=='k' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "NETWORK"

    // $ANTLR start "NTAX"
    public final void mNTAX() throws RecognitionException {
        try {
            int _type = NTAX;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:191:6: ( ( 'N' | 'n' ) ( 'T' | 't' ) ( 'A' | 'a' ) ( 'X' | 'x' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:191:8: ( 'N' | 'n' ) ( 'T' | 't' ) ( 'A' | 'a' ) ( 'X' | 'x' )
            {
            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='A'||input.LA(1)=='a' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='X'||input.LA(1)=='x' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "NTAX"

    // $ANTLR start "NCHAR"
    public final void mNCHAR() throws RecognitionException {
        try {
            int _type = NCHAR;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:193:7: ( ( 'N' | 'n' ) ( 'C' | 'c' ) ( 'H' | 'h' ) ( 'A' | 'a' ) ( 'R' | 'r' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:193:9: ( 'N' | 'n' ) ( 'C' | 'c' ) ( 'H' | 'h' ) ( 'A' | 'a' ) ( 'R' | 'r' )
            {
            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='C'||input.LA(1)=='c' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='H'||input.LA(1)=='h' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='A'||input.LA(1)=='a' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='R'||input.LA(1)=='r' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "NCHAR"

    // $ANTLR start "FORMAT"
    public final void mFORMAT() throws RecognitionException {
        try {
            int _type = FORMAT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:195:8: ( ( 'F' | 'f' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'M' | 'm' ) ( 'A' | 'a' ) ( 'T' | 't' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:195:10: ( 'F' | 'f' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'M' | 'm' ) ( 'A' | 'a' ) ( 'T' | 't' )
            {
            if ( input.LA(1)=='F'||input.LA(1)=='f' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='O'||input.LA(1)=='o' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='R'||input.LA(1)=='r' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='M'||input.LA(1)=='m' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='A'||input.LA(1)=='a' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "FORMAT"

    // $ANTLR start "DATATYPE"
    public final void mDATATYPE() throws RecognitionException {
        try {
            int _type = DATATYPE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:197:9: ( ( 'D' | 'd' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'Y' | 'y' ) ( 'P' | 'p' ) ( 'E' | 'e' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:197:11: ( 'D' | 'd' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'Y' | 'y' ) ( 'P' | 'p' ) ( 'E' | 'e' )
            {
            if ( input.LA(1)=='D'||input.LA(1)=='d' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='A'||input.LA(1)=='a' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='A'||input.LA(1)=='a' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='Y'||input.LA(1)=='y' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='P'||input.LA(1)=='p' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "DATATYPE"

    // $ANTLR start "SYMBOLS"
    public final void mSYMBOLS() throws RecognitionException {
        try {
            int _type = SYMBOLS;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:199:9: ( ( 'S' | 's' ) ( 'Y' | 'y' ) ( 'M' | 'm' ) ( 'B' | 'b' ) ( 'O' | 'o' ) ( 'L' | 'l' ) ( 'S' | 's' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:199:11: ( 'S' | 's' ) ( 'Y' | 'y' ) ( 'M' | 'm' ) ( 'B' | 'b' ) ( 'O' | 'o' ) ( 'L' | 'l' ) ( 'S' | 's' )
            {
            if ( input.LA(1)=='S'||input.LA(1)=='s' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='Y'||input.LA(1)=='y' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='M'||input.LA(1)=='m' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='B'||input.LA(1)=='b' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='O'||input.LA(1)=='o' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='L'||input.LA(1)=='l' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='S'||input.LA(1)=='s' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "SYMBOLS"

    // $ANTLR start "GAP"
    public final void mGAP() throws RecognitionException {
        try {
            int _type = GAP;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:201:5: ( ( 'G' | 'g' ) ( 'A' | 'a' ) ( 'P' | 'p' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:201:7: ( 'G' | 'g' ) ( 'A' | 'a' ) ( 'P' | 'p' )
            {
            if ( input.LA(1)=='G'||input.LA(1)=='g' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='A'||input.LA(1)=='a' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='P'||input.LA(1)=='p' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "GAP"

    // $ANTLR start "MATRIX"
    public final void mMATRIX() throws RecognitionException {
        try {
            int _type = MATRIX;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:203:8: ( ( 'M' | 'm' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'R' | 'r' ) ( 'I' | 'i' ) ( 'X' | 'x' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:203:10: ( 'M' | 'm' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'R' | 'r' ) ( 'I' | 'i' ) ( 'X' | 'x' )
            {
            if ( input.LA(1)=='M'||input.LA(1)=='m' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='A'||input.LA(1)=='a' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='R'||input.LA(1)=='r' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='I'||input.LA(1)=='i' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='X'||input.LA(1)=='x' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "MATRIX"

    // $ANTLR start "MISSING"
    public final void mMISSING() throws RecognitionException {
        try {
            int _type = MISSING;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:205:9: ( ( 'M' | 'm' ) ( 'I' | 'i' ) ( 'S' | 's' ) ( 'S' | 's' ) ( 'I' | 'i' ) ( 'N' | 'n' ) ( 'G' | 'g' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:205:11: ( 'M' | 'm' ) ( 'I' | 'i' ) ( 'S' | 's' ) ( 'S' | 's' ) ( 'I' | 'i' ) ( 'N' | 'n' ) ( 'G' | 'g' )
            {
            if ( input.LA(1)=='M'||input.LA(1)=='m' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='I'||input.LA(1)=='i' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='S'||input.LA(1)=='s' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='S'||input.LA(1)=='s' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='I'||input.LA(1)=='i' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='G'||input.LA(1)=='g' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "MISSING"

    // $ANTLR start "START"
    public final void mSTART() throws RecognitionException {
        try {
            int _type = START;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:208:2: ( '#' ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'X' | 'x' ) ( 'U' | 'u' ) ( 'S' | 's' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:208:4: '#' ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'X' | 'x' ) ( 'U' | 'u' ) ( 'S' | 's' )
            {
            match('#'); 

            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='X'||input.LA(1)=='x' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='U'||input.LA(1)=='u' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='S'||input.LA(1)=='s' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "START"

    // $ANTLR start "BEGIN"
    public final void mBEGIN() throws RecognitionException {
        try {
            int _type = BEGIN;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:210:8: ( ( 'B' | 'b' ) ( 'E' | 'e' ) ( 'G' | 'g' ) ( 'I' | 'i' ) ( 'N' | 'n' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:210:10: ( 'B' | 'b' ) ( 'E' | 'e' ) ( 'G' | 'g' ) ( 'I' | 'i' ) ( 'N' | 'n' )
            {
            if ( input.LA(1)=='B'||input.LA(1)=='b' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='G'||input.LA(1)=='g' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='I'||input.LA(1)=='i' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "BEGIN"

    // $ANTLR start "TRANSLATE"
    public final void mTRANSLATE() throws RecognitionException {
        try {
            int _type = TRANSLATE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:213:2: ( ( 'T' | 't' ) ( 'R' | 'r' ) ( 'A' | 'a' ) ( 'N' | 'n' ) ( 'S' | 's' ) ( 'L' | 'l' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'E' | 'e' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:213:4: ( 'T' | 't' ) ( 'R' | 'r' ) ( 'A' | 'a' ) ( 'N' | 'n' ) ( 'S' | 's' ) ( 'L' | 'l' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'E' | 'e' )
            {
            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='R'||input.LA(1)=='r' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='A'||input.LA(1)=='a' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='S'||input.LA(1)=='s' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='L'||input.LA(1)=='l' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='A'||input.LA(1)=='a' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "TRANSLATE"

    // $ANTLR start "NETWORKS"
    public final void mNETWORKS() throws RecognitionException {
        try {
            int _type = NETWORKS;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:216:2: ( ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' ) ( 'S' | 's' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:216:4: ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' ) ( 'S' | 's' )
            {
            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='W'||input.LA(1)=='w' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='O'||input.LA(1)=='o' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='R'||input.LA(1)=='r' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='K'||input.LA(1)=='k' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='S'||input.LA(1)=='s' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "NETWORKS"

    // $ANTLR start "TREES"
    public final void mTREES() throws RecognitionException {
        try {
            int _type = TREES;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:218:7: ( ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) ( 'S' | 's' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:218:9: ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) ( 'S' | 's' )
            {
            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='R'||input.LA(1)=='r' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='S'||input.LA(1)=='s' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "TREES"

    // $ANTLR start "END"
    public final void mEND() throws RecognitionException {
        try {
            int _type = END;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:220:5: ( ( 'E' | 'e' ) ( 'N' | 'n' ) ( 'D' | 'd' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:220:7: ( 'E' | 'e' ) ( 'N' | 'n' ) ( 'D' | 'd' )
            {
            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='D'||input.LA(1)=='d' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "END"

    // $ANTLR start "ID"
    public final void mID() throws RecognitionException {
        try {
            int _type = ID;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:222:4: ( ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' | '?' | '-' | '.' | '\\\\' | '/' )+ )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:222:6: ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' | '?' | '-' | '.' | '\\\\' | '/' )+
            {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:222:6: ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' | '?' | '-' | '.' | '\\\\' | '/' )+
            int cnt1=0;
            loop1:
            do {
                int alt1=2;
                int LA1_0 = input.LA(1);

                if ( ((LA1_0 >= '-' && LA1_0 <= '9')||LA1_0=='?'||(LA1_0 >= 'A' && LA1_0 <= 'Z')||LA1_0=='\\'||LA1_0=='_'||(LA1_0 >= 'a' && LA1_0 <= 'z')) ) {
                    alt1=1;
                }


                switch (alt1) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
            	    {
            	    if ( (input.LA(1) >= '-' && input.LA(1) <= '9')||input.LA(1)=='?'||(input.LA(1) >= 'A' && input.LA(1) <= 'Z')||input.LA(1)=='\\'||input.LA(1)=='_'||(input.LA(1) >= 'a' && input.LA(1) <= 'z') ) {
            	        input.consume();
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;
            	    }


            	    }
            	    break;

            	default :
            	    if ( cnt1 >= 1 ) break loop1;
                        EarlyExitException eee =
                            new EarlyExitException(1, input);
                        throw eee;
                }
                cnt1++;
            } while (true);


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "ID"

    // $ANTLR start "QUOTE"
    public final void mQUOTE() throws RecognitionException {
        try {
            int _type = QUOTE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:224:7: ( '\"' (~ ( '\"' ) )* '\"' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:224:9: '\"' (~ ( '\"' ) )* '\"'
            {
            match('\"'); 

            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:224:13: (~ ( '\"' ) )*
            loop2:
            do {
                int alt2=2;
                int LA2_0 = input.LA(1);

                if ( ((LA2_0 >= '\u0000' && LA2_0 <= '!')||(LA2_0 >= '#' && LA2_0 <= '\uFFFF')) ) {
                    alt2=1;
                }


                switch (alt2) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
            	    {
            	    if ( (input.LA(1) >= '\u0000' && input.LA(1) <= '!')||(input.LA(1) >= '#' && input.LA(1) <= '\uFFFF') ) {
            	        input.consume();
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;
            	    }


            	    }
            	    break;

            	default :
            	    break loop2;
                }
            } while (true);


            match('\"'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "QUOTE"

    // $ANTLR start "ID_SET"
    public final void mID_SET() throws RecognitionException {
        try {
            int _type = ID_SET;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:227:2: ( '{' ( ( WS )* ID ( WS )* ',' )* ( WS )* ID ( WS )* '}' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:227:4: '{' ( ( WS )* ID ( WS )* ',' )* ( WS )* ID ( WS )* '}'
            {
            match('{'); 

            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:227:8: ( ( WS )* ID ( WS )* ',' )*
            loop5:
            do {
                int alt5=2;
                alt5 = dfa5.predict(input);
                switch (alt5) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:227:9: ( WS )* ID ( WS )* ','
            	    {
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:227:9: ( WS )*
            	    loop3:
            	    do {
            	        int alt3=2;
            	        int LA3_0 = input.LA(1);

            	        if ( ((LA3_0 >= '\t' && LA3_0 <= '\n')||LA3_0=='\r'||LA3_0==' ') ) {
            	            alt3=1;
            	        }


            	        switch (alt3) {
            	    	case 1 :
            	    	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:227:9: WS
            	    	    {
            	    	    mWS(); 


            	    	    }
            	    	    break;

            	    	default :
            	    	    break loop3;
            	        }
            	    } while (true);


            	    mID(); 


            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:227:16: ( WS )*
            	    loop4:
            	    do {
            	        int alt4=2;
            	        int LA4_0 = input.LA(1);

            	        if ( ((LA4_0 >= '\t' && LA4_0 <= '\n')||LA4_0=='\r'||LA4_0==' ') ) {
            	            alt4=1;
            	        }


            	        switch (alt4) {
            	    	case 1 :
            	    	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:227:16: WS
            	    	    {
            	    	    mWS(); 


            	    	    }
            	    	    break;

            	    	default :
            	    	    break loop4;
            	        }
            	    } while (true);


            	    match(','); 

            	    }
            	    break;

            	default :
            	    break loop5;
                }
            } while (true);


            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:228:15: ( WS )*
            loop6:
            do {
                int alt6=2;
                int LA6_0 = input.LA(1);

                if ( ((LA6_0 >= '\t' && LA6_0 <= '\n')||LA6_0=='\r'||LA6_0==' ') ) {
                    alt6=1;
                }


                switch (alt6) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:228:15: WS
            	    {
            	    mWS(); 


            	    }
            	    break;

            	default :
            	    break loop6;
                }
            } while (true);


            mID(); 


            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:228:22: ( WS )*
            loop7:
            do {
                int alt7=2;
                int LA7_0 = input.LA(1);

                if ( ((LA7_0 >= '\t' && LA7_0 <= '\n')||LA7_0=='\r'||LA7_0==' ') ) {
                    alt7=1;
                }


                switch (alt7) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:228:22: WS
            	    {
            	    mWS(); 


            	    }
            	    break;

            	default :
            	    break loop7;
                }
            } while (true);


            match('}'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "ID_SET"

    // $ANTLR start "TAXON_SET_LIST"
    public final void mTAXON_SET_LIST() throws RecognitionException {
        try {
            int _type = TAXON_SET_LIST;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:231:2: ( '(' ( ( WS )* ID_SET ( ( WS )* ',' )? )* ( ( WS )* ID_SET ( WS )* ) ')' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:231:4: '(' ( ( WS )* ID_SET ( ( WS )* ',' )? )* ( ( WS )* ID_SET ( WS )* ) ')'
            {
            match('('); 

            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:231:8: ( ( WS )* ID_SET ( ( WS )* ',' )? )*
            loop11:
            do {
                int alt11=2;
                alt11 = dfa11.predict(input);
                switch (alt11) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:231:9: ( WS )* ID_SET ( ( WS )* ',' )?
            	    {
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:231:9: ( WS )*
            	    loop8:
            	    do {
            	        int alt8=2;
            	        int LA8_0 = input.LA(1);

            	        if ( ((LA8_0 >= '\t' && LA8_0 <= '\n')||LA8_0=='\r'||LA8_0==' ') ) {
            	            alt8=1;
            	        }


            	        switch (alt8) {
            	    	case 1 :
            	    	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:231:9: WS
            	    	    {
            	    	    mWS(); 


            	    	    }
            	    	    break;

            	    	default :
            	    	    break loop8;
            	        }
            	    } while (true);


            	    mID_SET(); 


            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:231:20: ( ( WS )* ',' )?
            	    int alt10=2;
            	    alt10 = dfa10.predict(input);
            	    switch (alt10) {
            	        case 1 :
            	            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:231:21: ( WS )* ','
            	            {
            	            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:231:21: ( WS )*
            	            loop9:
            	            do {
            	                int alt9=2;
            	                int LA9_0 = input.LA(1);

            	                if ( ((LA9_0 >= '\t' && LA9_0 <= '\n')||LA9_0=='\r'||LA9_0==' ') ) {
            	                    alt9=1;
            	                }


            	                switch (alt9) {
            	            	case 1 :
            	            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:231:21: WS
            	            	    {
            	            	    mWS(); 


            	            	    }
            	            	    break;

            	            	default :
            	            	    break loop9;
            	                }
            	            } while (true);


            	            match(','); 

            	            }
            	            break;

            	    }


            	    }
            	    break;

            	default :
            	    break loop11;
                }
            } while (true);


            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:232:14: ( ( WS )* ID_SET ( WS )* )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:232:15: ( WS )* ID_SET ( WS )*
            {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:232:15: ( WS )*
            loop12:
            do {
                int alt12=2;
                int LA12_0 = input.LA(1);

                if ( ((LA12_0 >= '\t' && LA12_0 <= '\n')||LA12_0=='\r'||LA12_0==' ') ) {
                    alt12=1;
                }


                switch (alt12) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:232:15: WS
            	    {
            	    mWS(); 


            	    }
            	    break;

            	default :
            	    break loop12;
                }
            } while (true);


            mID_SET(); 


            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:232:26: ( WS )*
            loop13:
            do {
                int alt13=2;
                int LA13_0 = input.LA(1);

                if ( ((LA13_0 >= '\t' && LA13_0 <= '\n')||LA13_0=='\r'||LA13_0==' ') ) {
                    alt13=1;
                }


                switch (alt13) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:232:26: WS
            	    {
            	    mWS(); 


            	    }
            	    break;

            	default :
            	    break loop13;
                }
            } while (true);


            }


            match(')'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "TAXON_SET_LIST"

    // $ANTLR start "RN_LS_NONCOMMENT"
    public final void mRN_LS_NONCOMMENT() throws RecognitionException {
        try {
            int _type = RN_LS_NONCOMMENT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:235:2: ( '[' '&' (~ ( ']' ) )* ']' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:235:4: '[' '&' (~ ( ']' ) )* ']'
            {
            match('['); 

            match('&'); 

            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:235:12: (~ ( ']' ) )*
            loop14:
            do {
                int alt14=2;
                int LA14_0 = input.LA(1);

                if ( ((LA14_0 >= '\u0000' && LA14_0 <= '\\')||(LA14_0 >= '^' && LA14_0 <= '\uFFFF')) ) {
                    alt14=1;
                }


                switch (alt14) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
            	    {
            	    if ( (input.LA(1) >= '\u0000' && input.LA(1) <= '\\')||(input.LA(1) >= '^' && input.LA(1) <= '\uFFFF') ) {
            	        input.consume();
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;
            	    }


            	    }
            	    break;

            	default :
            	    break loop14;
                }
            } while (true);


            match(']'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "RN_LS_NONCOMMENT"

    // $ANTLR start "NESTED_ML_COMMENT"
    public final void mNESTED_ML_COMMENT() throws RecognitionException {
        try {
            int _type = NESTED_ML_COMMENT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:238:3: ( '[' (~ ( '[' | ']' ) )* ']' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:238:6: '[' (~ ( '[' | ']' ) )* ']'
            {
            match('['); 

            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:238:10: (~ ( '[' | ']' ) )*
            loop15:
            do {
                int alt15=2;
                int LA15_0 = input.LA(1);

                if ( ((LA15_0 >= '\u0000' && LA15_0 <= 'Z')||LA15_0=='\\'||(LA15_0 >= '^' && LA15_0 <= '\uFFFF')) ) {
                    alt15=1;
                }


                switch (alt15) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
            	    {
            	    if ( (input.LA(1) >= '\u0000' && input.LA(1) <= 'Z')||input.LA(1)=='\\'||(input.LA(1) >= '^' && input.LA(1) <= '\uFFFF') ) {
            	        input.consume();
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;
            	    }


            	    }
            	    break;

            	default :
            	    break loop15;
                }
            } while (true);


            match(']'); 

            _channel=HIDDEN;

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "NESTED_ML_COMMENT"

    // $ANTLR start "WS"
    public final void mWS() throws RecognitionException {
        try {
            int _type = WS;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:240:5: ( ( ' ' | '\\t' | '\\r' | '\\n' ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:240:9: ( ' ' | '\\t' | '\\r' | '\\n' )
            {
            if ( (input.LA(1) >= '\t' && input.LA(1) <= '\n')||input.LA(1)=='\r'||input.LA(1)==' ' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            _channel=HIDDEN;

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "WS"

    // $ANTLR start "ELSE"
    public final void mELSE() throws RecognitionException {
        try {
            int _type = ELSE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:249:6: ( . )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:249:8: .
            {
            matchAny(); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "ELSE"

    public void mTokens() throws RecognitionException {
        // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:8: ( T__34 | T__35 | T__36 | T__37 | T__38 | T__39 | T__40 | T__41 | DEFAULT_INDICATOR | PHYLONET | DATA | MORPHDATA | DIMENSIONS | TREE | UTREE | NETWORK | NTAX | NCHAR | FORMAT | DATATYPE | SYMBOLS | GAP | MATRIX | MISSING | START | BEGIN | TRANSLATE | NETWORKS | TREES | END | ID | QUOTE | ID_SET | TAXON_SET_LIST | RN_LS_NONCOMMENT | NESTED_ML_COMMENT | WS | ELSE )
        int alt16=38;
        alt16 = dfa16.predict(input);
        switch (alt16) {
            case 1 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:10: T__34
                {
                mT__34(); 


                }
                break;
            case 2 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:16: T__35
                {
                mT__35(); 


                }
                break;
            case 3 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:22: T__36
                {
                mT__36(); 


                }
                break;
            case 4 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:28: T__37
                {
                mT__37(); 


                }
                break;
            case 5 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:34: T__38
                {
                mT__38(); 


                }
                break;
            case 6 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:40: T__39
                {
                mT__39(); 


                }
                break;
            case 7 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:46: T__40
                {
                mT__40(); 


                }
                break;
            case 8 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:52: T__41
                {
                mT__41(); 


                }
                break;
            case 9 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:58: DEFAULT_INDICATOR
                {
                mDEFAULT_INDICATOR(); 


                }
                break;
            case 10 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:76: PHYLONET
                {
                mPHYLONET(); 


                }
                break;
            case 11 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:85: DATA
                {
                mDATA(); 


                }
                break;
            case 12 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:90: MORPHDATA
                {
                mMORPHDATA(); 


                }
                break;
            case 13 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:100: DIMENSIONS
                {
                mDIMENSIONS(); 


                }
                break;
            case 14 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:111: TREE
                {
                mTREE(); 


                }
                break;
            case 15 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:116: UTREE
                {
                mUTREE(); 


                }
                break;
            case 16 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:122: NETWORK
                {
                mNETWORK(); 


                }
                break;
            case 17 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:130: NTAX
                {
                mNTAX(); 


                }
                break;
            case 18 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:135: NCHAR
                {
                mNCHAR(); 


                }
                break;
            case 19 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:141: FORMAT
                {
                mFORMAT(); 


                }
                break;
            case 20 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:148: DATATYPE
                {
                mDATATYPE(); 


                }
                break;
            case 21 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:157: SYMBOLS
                {
                mSYMBOLS(); 


                }
                break;
            case 22 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:165: GAP
                {
                mGAP(); 


                }
                break;
            case 23 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:169: MATRIX
                {
                mMATRIX(); 


                }
                break;
            case 24 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:176: MISSING
                {
                mMISSING(); 


                }
                break;
            case 25 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:184: START
                {
                mSTART(); 


                }
                break;
            case 26 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:190: BEGIN
                {
                mBEGIN(); 


                }
                break;
            case 27 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:196: TRANSLATE
                {
                mTRANSLATE(); 


                }
                break;
            case 28 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:206: NETWORKS
                {
                mNETWORKS(); 


                }
                break;
            case 29 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:215: TREES
                {
                mTREES(); 


                }
                break;
            case 30 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:221: END
                {
                mEND(); 


                }
                break;
            case 31 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:225: ID
                {
                mID(); 


                }
                break;
            case 32 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:228: QUOTE
                {
                mQUOTE(); 


                }
                break;
            case 33 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:234: ID_SET
                {
                mID_SET(); 


                }
                break;
            case 34 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:241: TAXON_SET_LIST
                {
                mTAXON_SET_LIST(); 


                }
                break;
            case 35 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:256: RN_LS_NONCOMMENT
                {
                mRN_LS_NONCOMMENT(); 


                }
                break;
            case 36 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:273: NESTED_ML_COMMENT
                {
                mNESTED_ML_COMMENT(); 


                }
                break;
            case 37 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:291: WS
                {
                mWS(); 


                }
                break;
            case 38 :
                // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:294: ELSE
                {
                mELSE(); 


                }
                break;

        }

    }


    protected DFA5 dfa5 = new DFA5(this);
    protected DFA11 dfa11 = new DFA11(this);
    protected DFA10 dfa10 = new DFA10(this);
    protected DFA16 dfa16 = new DFA16(this);
    static final String DFA5_eotS =
        "\6\uffff";
    static final String DFA5_eofS =
        "\6\uffff";
    static final String DFA5_minS =
        "\4\11\2\uffff";
    static final String DFA5_maxS =
        "\2\172\2\175\2\uffff";
    static final String DFA5_acceptS =
        "\4\uffff\1\2\1\1";
    static final String DFA5_specialS =
        "\6\uffff}>";
    static final String[] DFA5_transitionS = {
            "\2\1\2\uffff\1\1\22\uffff\1\1\14\uffff\15\2\5\uffff\1\2\1\uffff"+
            "\32\2\1\uffff\1\2\2\uffff\1\2\1\uffff\32\2",
            "\2\1\2\uffff\1\1\22\uffff\1\1\14\uffff\15\2\5\uffff\1\2\1\uffff"+
            "\32\2\1\uffff\1\2\2\uffff\1\2\1\uffff\32\2",
            "\2\3\2\uffff\1\3\22\uffff\1\3\13\uffff\1\5\15\2\5\uffff\1\2"+
            "\1\uffff\32\2\1\uffff\1\2\2\uffff\1\2\1\uffff\32\2\2\uffff\1"+
            "\4",
            "\2\3\2\uffff\1\3\22\uffff\1\3\13\uffff\1\5\120\uffff\1\4",
            "",
            ""
    };

    static final short[] DFA5_eot = DFA.unpackEncodedString(DFA5_eotS);
    static final short[] DFA5_eof = DFA.unpackEncodedString(DFA5_eofS);
    static final char[] DFA5_min = DFA.unpackEncodedStringToUnsignedChars(DFA5_minS);
    static final char[] DFA5_max = DFA.unpackEncodedStringToUnsignedChars(DFA5_maxS);
    static final short[] DFA5_accept = DFA.unpackEncodedString(DFA5_acceptS);
    static final short[] DFA5_special = DFA.unpackEncodedString(DFA5_specialS);
    static final short[][] DFA5_transition;

    static {
        int numStates = DFA5_transitionS.length;
        DFA5_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA5_transition[i] = DFA.unpackEncodedString(DFA5_transitionS[i]);
        }
    }

    class DFA5 extends DFA {

        public DFA5(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 5;
            this.eot = DFA5_eot;
            this.eof = DFA5_eof;
            this.min = DFA5_min;
            this.max = DFA5_max;
            this.accept = DFA5_accept;
            this.special = DFA5_special;
            this.transition = DFA5_transition;
        }
        public String getDescription() {
            return "()* loopback of 227:8: ( ( WS )* ID ( WS )* ',' )*";
        }
    }
    static final String DFA11_eotS =
        "\13\uffff";
    static final String DFA11_eofS =
        "\13\uffff";
    static final String DFA11_minS =
        "\11\11\2\uffff";
    static final String DFA11_maxS =
        "\2\173\2\172\2\175\1\172\2\173\2\uffff";
    static final String DFA11_acceptS =
        "\11\uffff\1\2\1\1";
    static final String DFA11_specialS =
        "\13\uffff}>";
    static final String[] DFA11_transitionS = {
            "\2\1\2\uffff\1\1\22\uffff\1\1\132\uffff\1\2",
            "\2\1\2\uffff\1\1\22\uffff\1\1\132\uffff\1\2",
            "\2\3\2\uffff\1\3\22\uffff\1\3\14\uffff\15\4\5\uffff\1\4\1\uffff"+
            "\32\4\1\uffff\1\4\2\uffff\1\4\1\uffff\32\4",
            "\2\3\2\uffff\1\3\22\uffff\1\3\14\uffff\15\4\5\uffff\1\4\1\uffff"+
            "\32\4\1\uffff\1\4\2\uffff\1\4\1\uffff\32\4",
            "\2\5\2\uffff\1\5\22\uffff\1\5\13\uffff\1\6\15\4\5\uffff\1\4"+
            "\1\uffff\32\4\1\uffff\1\4\2\uffff\1\4\1\uffff\32\4\2\uffff\1"+
            "\7",
            "\2\5\2\uffff\1\5\22\uffff\1\5\13\uffff\1\6\120\uffff\1\7",
            "\2\3\2\uffff\1\3\22\uffff\1\3\14\uffff\15\4\5\uffff\1\4\1\uffff"+
            "\32\4\1\uffff\1\4\2\uffff\1\4\1\uffff\32\4",
            "\2\10\2\uffff\1\10\22\uffff\1\10\10\uffff\1\11\2\uffff\1\12"+
            "\116\uffff\1\12",
            "\2\10\2\uffff\1\10\22\uffff\1\10\10\uffff\1\11\2\uffff\1\12"+
            "\116\uffff\1\12",
            "",
            ""
    };

    static final short[] DFA11_eot = DFA.unpackEncodedString(DFA11_eotS);
    static final short[] DFA11_eof = DFA.unpackEncodedString(DFA11_eofS);
    static final char[] DFA11_min = DFA.unpackEncodedStringToUnsignedChars(DFA11_minS);
    static final char[] DFA11_max = DFA.unpackEncodedStringToUnsignedChars(DFA11_maxS);
    static final short[] DFA11_accept = DFA.unpackEncodedString(DFA11_acceptS);
    static final short[] DFA11_special = DFA.unpackEncodedString(DFA11_specialS);
    static final short[][] DFA11_transition;

    static {
        int numStates = DFA11_transitionS.length;
        DFA11_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA11_transition[i] = DFA.unpackEncodedString(DFA11_transitionS[i]);
        }
    }

    class DFA11 extends DFA {

        public DFA11(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 11;
            this.eot = DFA11_eot;
            this.eof = DFA11_eof;
            this.min = DFA11_min;
            this.max = DFA11_max;
            this.accept = DFA11_accept;
            this.special = DFA11_special;
            this.transition = DFA11_transition;
        }
        public String getDescription() {
            return "()* loopback of 231:8: ( ( WS )* ID_SET ( ( WS )* ',' )? )*";
        }
    }
    static final String DFA10_eotS =
        "\4\uffff";
    static final String DFA10_eofS =
        "\4\uffff";
    static final String DFA10_minS =
        "\2\11\2\uffff";
    static final String DFA10_maxS =
        "\2\173\2\uffff";
    static final String DFA10_acceptS =
        "\2\uffff\1\1\1\2";
    static final String DFA10_specialS =
        "\4\uffff}>";
    static final String[] DFA10_transitionS = {
            "\2\1\2\uffff\1\1\22\uffff\1\1\13\uffff\1\2\116\uffff\1\3",
            "\2\1\2\uffff\1\1\22\uffff\1\1\13\uffff\1\2\116\uffff\1\3",
            "",
            ""
    };

    static final short[] DFA10_eot = DFA.unpackEncodedString(DFA10_eotS);
    static final short[] DFA10_eof = DFA.unpackEncodedString(DFA10_eofS);
    static final char[] DFA10_min = DFA.unpackEncodedStringToUnsignedChars(DFA10_minS);
    static final char[] DFA10_max = DFA.unpackEncodedStringToUnsignedChars(DFA10_maxS);
    static final short[] DFA10_accept = DFA.unpackEncodedString(DFA10_acceptS);
    static final short[] DFA10_special = DFA.unpackEncodedString(DFA10_specialS);
    static final short[][] DFA10_transition;

    static {
        int numStates = DFA10_transitionS.length;
        DFA10_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA10_transition[i] = DFA.unpackEncodedString(DFA10_transitionS[i]);
        }
    }

    class DFA10 extends DFA {

        public DFA10(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 10;
            this.eot = DFA10_eot;
            this.eof = DFA10_eof;
            this.min = DFA10_min;
            this.max = DFA10_max;
            this.accept = DFA10_accept;
            this.special = DFA10_special;
            this.transition = DFA10_transition;
        }
        public String getDescription() {
            return "231:20: ( ( WS )* ',' )?";
        }
    }
    static final String DFA16_eotS =
        "\1\uffff\1\34\10\uffff\11\47\1\33\2\47\1\uffff\3\33\14\uffff\1\47"+
        "\1\uffff\15\47\1\uffff\2\47\5\uffff\16\47\1\137\1\47\1\141\3\uffff"+
        "\1\47\1\144\4\47\1\152\3\47\1\156\3\47\1\uffff\1\47\1\uffff\2\47"+
        "\1\uffff\4\47\1\171\1\uffff\1\47\1\173\1\47\1\uffff\1\175\2\47\1"+
        "\u0080\4\47\1\u0085\1\47\1\uffff\1\47\1\uffff\1\47\1\uffff\1\u0089"+
        "\1\47\1\uffff\4\47\1\uffff\1\u008f\1\47\1\u0092\1\uffff\1\u0093"+
        "\1\u0094\1\u0095\2\47\1\uffff\1\47\1\u0099\4\uffff\1\47\1\u009b"+
        "\1\u009c\1\uffff\1\u009d\3\uffff";
    static final String DFA16_eofS =
        "\u009e\uffff";
    static final String DFA16_minS =
        "\1\0\1\11\10\uffff\1\110\2\101\1\122\1\124\1\103\1\117\1\131\1\101"+
        "\1\116\1\105\1\116\1\uffff\1\0\1\11\1\0\14\uffff\1\131\1\uffff\1"+
        "\124\1\115\1\122\1\124\1\123\1\101\1\122\1\124\1\101\1\110\1\122"+
        "\1\115\1\120\1\uffff\1\107\1\104\2\uffff\1\0\2\uffff\1\114\1\101"+
        "\1\105\1\120\1\122\1\123\1\105\1\116\1\105\1\127\1\130\1\101\1\115"+
        "\1\102\1\55\1\111\1\55\1\0\2\uffff\1\117\1\55\1\116\1\110\2\111"+
        "\1\55\1\123\1\105\1\117\1\55\1\122\1\101\1\117\1\uffff\1\116\1\uffff"+
        "\1\116\1\131\1\uffff\1\123\1\104\1\130\1\116\1\55\1\uffff\1\114"+
        "\1\55\1\122\1\uffff\1\55\1\124\1\114\1\55\1\105\1\120\1\111\1\101"+
        "\1\55\1\107\1\uffff\1\101\1\uffff\1\113\1\uffff\1\55\1\123\1\uffff"+
        "\1\124\1\105\1\117\1\124\1\uffff\1\55\1\124\1\55\1\uffff\3\55\1"+
        "\116\1\101\1\uffff\1\105\1\55\4\uffff\1\123\2\55\1\uffff\1\55\3"+
        "\uffff";
    static final String DFA16_maxS =
        "\1\uffff\1\173\10\uffff\1\150\1\151\1\157\1\162\2\164\1\157\1\171"+
        "\1\141\1\156\1\145\1\156\1\uffff\1\uffff\1\172\1\uffff\14\uffff"+
        "\1\171\1\uffff\1\164\1\155\1\162\1\164\1\163\1\145\1\162\1\164\1"+
        "\141\1\150\1\162\1\155\1\160\1\uffff\1\147\1\144\2\uffff\1\uffff"+
        "\2\uffff\1\154\1\141\1\145\1\160\1\162\1\163\1\145\1\156\1\145\1"+
        "\167\1\170\1\141\1\155\1\142\1\172\1\151\1\172\1\uffff\2\uffff\1"+
        "\157\1\172\1\156\1\150\2\151\1\172\1\163\1\145\1\157\1\172\1\162"+
        "\1\141\1\157\1\uffff\1\156\1\uffff\1\156\1\171\1\uffff\1\163\1\144"+
        "\1\170\1\156\1\172\1\uffff\1\154\1\172\1\162\1\uffff\1\172\1\164"+
        "\1\154\1\172\1\145\1\160\1\151\1\141\1\172\1\147\1\uffff\1\141\1"+
        "\uffff\1\153\1\uffff\1\172\1\163\1\uffff\1\164\1\145\1\157\1\164"+
        "\1\uffff\1\172\1\164\1\172\1\uffff\3\172\1\156\1\141\1\uffff\1\145"+
        "\1\172\4\uffff\1\163\2\172\1\uffff\1\172\3\uffff";
    static final String DFA16_acceptS =
        "\2\uffff\1\2\1\3\1\4\1\5\1\6\1\7\1\10\1\11\14\uffff\1\37\3\uffff"+
        "\1\45\1\46\1\1\1\42\1\2\1\3\1\4\1\5\1\6\1\7\1\10\1\11\1\uffff\1"+
        "\37\15\uffff\1\31\2\uffff\1\40\1\41\1\uffff\1\44\1\45\22\uffff\2"+
        "\43\16\uffff\1\26\1\uffff\1\36\2\uffff\1\13\5\uffff\1\16\3\uffff"+
        "\1\21\12\uffff\1\35\1\uffff\1\17\1\uffff\1\22\2\uffff\1\32\4\uffff"+
        "\1\27\3\uffff\1\23\5\uffff\1\30\2\uffff\1\20\1\25\1\12\1\24\3\uffff"+
        "\1\34\1\uffff\1\14\1\33\1\15";
    static final String DFA16_specialS =
        "\1\4\26\uffff\1\2\1\uffff\1\0\40\uffff\1\1\23\uffff\1\3\117\uffff}>";
    static final String[] DFA16_transitionS = {
            "\11\33\2\32\2\33\1\32\22\33\1\32\1\33\1\27\1\23\4\33\1\1\1\2"+
            "\1\11\1\33\1\3\15\26\1\4\1\5\1\6\1\7\1\10\1\26\1\33\1\26\1\24"+
            "\1\26\1\13\1\25\1\20\1\22\5\26\1\14\1\17\1\26\1\12\2\26\1\21"+
            "\1\15\1\16\5\26\1\31\1\26\2\33\1\26\1\33\1\26\1\24\1\26\1\13"+
            "\1\25\1\20\1\22\5\26\1\14\1\17\1\26\1\12\2\26\1\21\1\15\1\16"+
            "\5\26\1\30\uff84\33",
            "\2\35\2\uffff\1\35\22\uffff\1\35\132\uffff\1\35",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "\1\46\37\uffff\1\46",
            "\1\50\7\uffff\1\51\27\uffff\1\50\7\uffff\1\51",
            "\1\53\7\uffff\1\54\5\uffff\1\52\21\uffff\1\53\7\uffff\1\54"+
            "\5\uffff\1\52",
            "\1\55\37\uffff\1\55",
            "\1\56\37\uffff\1\56",
            "\1\61\1\uffff\1\57\16\uffff\1\60\16\uffff\1\61\1\uffff\1\57"+
            "\16\uffff\1\60",
            "\1\62\37\uffff\1\62",
            "\1\63\37\uffff\1\63",
            "\1\64\37\uffff\1\64",
            "\1\65\37\uffff\1\65",
            "\1\66\37\uffff\1\66",
            "\1\67\37\uffff\1\67",
            "",
            "\0\70",
            "\2\71\2\uffff\1\71\22\uffff\1\71\14\uffff\15\71\5\uffff\1\71"+
            "\1\uffff\32\71\1\uffff\1\71\2\uffff\1\71\1\uffff\32\71",
            "\46\73\1\72\64\73\1\uffff\uffa4\73",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "\1\75\37\uffff\1\75",
            "",
            "\1\76\37\uffff\1\76",
            "\1\77\37\uffff\1\77",
            "\1\100\37\uffff\1\100",
            "\1\101\37\uffff\1\101",
            "\1\102\37\uffff\1\102",
            "\1\104\3\uffff\1\103\33\uffff\1\104\3\uffff\1\103",
            "\1\105\37\uffff\1\105",
            "\1\106\37\uffff\1\106",
            "\1\107\37\uffff\1\107",
            "\1\110\37\uffff\1\110",
            "\1\111\37\uffff\1\111",
            "\1\112\37\uffff\1\112",
            "\1\113\37\uffff\1\113",
            "",
            "\1\114\37\uffff\1\114",
            "\1\115\37\uffff\1\115",
            "",
            "",
            "\133\116\1\120\1\116\1\117\uffa2\116",
            "",
            "",
            "\1\121\37\uffff\1\121",
            "\1\122\37\uffff\1\122",
            "\1\123\37\uffff\1\123",
            "\1\124\37\uffff\1\124",
            "\1\125\37\uffff\1\125",
            "\1\126\37\uffff\1\126",
            "\1\127\37\uffff\1\127",
            "\1\130\37\uffff\1\130",
            "\1\131\37\uffff\1\131",
            "\1\132\37\uffff\1\132",
            "\1\133\37\uffff\1\133",
            "\1\134\37\uffff\1\134",
            "\1\135\37\uffff\1\135",
            "\1\136\37\uffff\1\136",
            "\15\47\5\uffff\1\47\1\uffff\32\47\1\uffff\1\47\2\uffff\1\47"+
            "\1\uffff\32\47",
            "\1\140\37\uffff\1\140",
            "\15\47\5\uffff\1\47\1\uffff\32\47\1\uffff\1\47\2\uffff\1\47"+
            "\1\uffff\32\47",
            "\133\116\1\120\1\116\1\117\uffa2\116",
            "",
            "",
            "\1\142\37\uffff\1\142",
            "\15\47\5\uffff\1\47\1\uffff\23\47\1\143\6\47\1\uffff\1\47\2"+
            "\uffff\1\47\1\uffff\23\47\1\143\6\47",
            "\1\145\37\uffff\1\145",
            "\1\146\37\uffff\1\146",
            "\1\147\37\uffff\1\147",
            "\1\150\37\uffff\1\150",
            "\15\47\5\uffff\1\47\1\uffff\22\47\1\151\7\47\1\uffff\1\47\2"+
            "\uffff\1\47\1\uffff\22\47\1\151\7\47",
            "\1\153\37\uffff\1\153",
            "\1\154\37\uffff\1\154",
            "\1\155\37\uffff\1\155",
            "\15\47\5\uffff\1\47\1\uffff\32\47\1\uffff\1\47\2\uffff\1\47"+
            "\1\uffff\32\47",
            "\1\157\37\uffff\1\157",
            "\1\160\37\uffff\1\160",
            "\1\161\37\uffff\1\161",
            "",
            "\1\162\37\uffff\1\162",
            "",
            "\1\163\37\uffff\1\163",
            "\1\164\37\uffff\1\164",
            "",
            "\1\165\37\uffff\1\165",
            "\1\166\37\uffff\1\166",
            "\1\167\37\uffff\1\167",
            "\1\170\37\uffff\1\170",
            "\15\47\5\uffff\1\47\1\uffff\32\47\1\uffff\1\47\2\uffff\1\47"+
            "\1\uffff\32\47",
            "",
            "\1\172\37\uffff\1\172",
            "\15\47\5\uffff\1\47\1\uffff\32\47\1\uffff\1\47\2\uffff\1\47"+
            "\1\uffff\32\47",
            "\1\174\37\uffff\1\174",
            "",
            "\15\47\5\uffff\1\47\1\uffff\32\47\1\uffff\1\47\2\uffff\1\47"+
            "\1\uffff\32\47",
            "\1\176\37\uffff\1\176",
            "\1\177\37\uffff\1\177",
            "\15\47\5\uffff\1\47\1\uffff\32\47\1\uffff\1\47\2\uffff\1\47"+
            "\1\uffff\32\47",
            "\1\u0081\37\uffff\1\u0081",
            "\1\u0082\37\uffff\1\u0082",
            "\1\u0083\37\uffff\1\u0083",
            "\1\u0084\37\uffff\1\u0084",
            "\15\47\5\uffff\1\47\1\uffff\32\47\1\uffff\1\47\2\uffff\1\47"+
            "\1\uffff\32\47",
            "\1\u0086\37\uffff\1\u0086",
            "",
            "\1\u0087\37\uffff\1\u0087",
            "",
            "\1\u0088\37\uffff\1\u0088",
            "",
            "\15\47\5\uffff\1\47\1\uffff\32\47\1\uffff\1\47\2\uffff\1\47"+
            "\1\uffff\32\47",
            "\1\u008a\37\uffff\1\u008a",
            "",
            "\1\u008b\37\uffff\1\u008b",
            "\1\u008c\37\uffff\1\u008c",
            "\1\u008d\37\uffff\1\u008d",
            "\1\u008e\37\uffff\1\u008e",
            "",
            "\15\47\5\uffff\1\47\1\uffff\32\47\1\uffff\1\47\2\uffff\1\47"+
            "\1\uffff\32\47",
            "\1\u0090\37\uffff\1\u0090",
            "\15\47\5\uffff\1\47\1\uffff\22\47\1\u0091\7\47\1\uffff\1\47"+
            "\2\uffff\1\47\1\uffff\22\47\1\u0091\7\47",
            "",
            "\15\47\5\uffff\1\47\1\uffff\32\47\1\uffff\1\47\2\uffff\1\47"+
            "\1\uffff\32\47",
            "\15\47\5\uffff\1\47\1\uffff\32\47\1\uffff\1\47\2\uffff\1\47"+
            "\1\uffff\32\47",
            "\15\47\5\uffff\1\47\1\uffff\32\47\1\uffff\1\47\2\uffff\1\47"+
            "\1\uffff\32\47",
            "\1\u0096\37\uffff\1\u0096",
            "\1\u0097\37\uffff\1\u0097",
            "",
            "\1\u0098\37\uffff\1\u0098",
            "\15\47\5\uffff\1\47\1\uffff\32\47\1\uffff\1\47\2\uffff\1\47"+
            "\1\uffff\32\47",
            "",
            "",
            "",
            "",
            "\1\u009a\37\uffff\1\u009a",
            "\15\47\5\uffff\1\47\1\uffff\32\47\1\uffff\1\47\2\uffff\1\47"+
            "\1\uffff\32\47",
            "\15\47\5\uffff\1\47\1\uffff\32\47\1\uffff\1\47\2\uffff\1\47"+
            "\1\uffff\32\47",
            "",
            "\15\47\5\uffff\1\47\1\uffff\32\47\1\uffff\1\47\2\uffff\1\47"+
            "\1\uffff\32\47",
            "",
            "",
            ""
    };

    static final short[] DFA16_eot = DFA.unpackEncodedString(DFA16_eotS);
    static final short[] DFA16_eof = DFA.unpackEncodedString(DFA16_eofS);
    static final char[] DFA16_min = DFA.unpackEncodedStringToUnsignedChars(DFA16_minS);
    static final char[] DFA16_max = DFA.unpackEncodedStringToUnsignedChars(DFA16_maxS);
    static final short[] DFA16_accept = DFA.unpackEncodedString(DFA16_acceptS);
    static final short[] DFA16_special = DFA.unpackEncodedString(DFA16_specialS);
    static final short[][] DFA16_transition;

    static {
        int numStates = DFA16_transitionS.length;
        DFA16_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA16_transition[i] = DFA.unpackEncodedString(DFA16_transitionS[i]);
        }
    }

    class DFA16 extends DFA {

        public DFA16(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 16;
            this.eot = DFA16_eot;
            this.eof = DFA16_eof;
            this.min = DFA16_min;
            this.max = DFA16_max;
            this.accept = DFA16_accept;
            this.special = DFA16_special;
            this.transition = DFA16_transition;
        }
        public String getDescription() {
            return "1:1: Tokens : ( T__34 | T__35 | T__36 | T__37 | T__38 | T__39 | T__40 | T__41 | DEFAULT_INDICATOR | PHYLONET | DATA | MORPHDATA | DIMENSIONS | TREE | UTREE | NETWORK | NTAX | NCHAR | FORMAT | DATATYPE | SYMBOLS | GAP | MATRIX | MISSING | START | BEGIN | TRANSLATE | NETWORKS | TREES | END | ID | QUOTE | ID_SET | TAXON_SET_LIST | RN_LS_NONCOMMENT | NESTED_ML_COMMENT | WS | ELSE );";
        }
        public int specialStateTransition(int s, IntStream _input) throws NoViableAltException {
            IntStream input = _input;
        	int _s = s;
            switch ( s ) {
                    case 0 : 
                        int LA16_25 = input.LA(1);

                        s = -1;
                        if ( (LA16_25=='&') ) {s = 58;}

                        else if ( ((LA16_25 >= '\u0000' && LA16_25 <= '%')||(LA16_25 >= '\'' && LA16_25 <= 'Z')||(LA16_25 >= '\\' && LA16_25 <= '\uFFFF')) ) {s = 59;}

                        else s = 27;

                        if ( s>=0 ) return s;
                        break;

                    case 1 : 
                        int LA16_58 = input.LA(1);

                        s = -1;
                        if ( ((LA16_58 >= '\u0000' && LA16_58 <= 'Z')||LA16_58=='\\'||(LA16_58 >= '^' && LA16_58 <= '\uFFFF')) ) {s = 78;}

                        else if ( (LA16_58==']') ) {s = 79;}

                        else if ( (LA16_58=='[') ) {s = 80;}

                        if ( s>=0 ) return s;
                        break;

                    case 2 : 
                        int LA16_23 = input.LA(1);

                        s = -1;
                        if ( ((LA16_23 >= '\u0000' && LA16_23 <= '\uFFFF')) ) {s = 56;}

                        else s = 27;

                        if ( s>=0 ) return s;
                        break;

                    case 3 : 
                        int LA16_78 = input.LA(1);

                        s = -1;
                        if ( (LA16_78==']') ) {s = 79;}

                        else if ( ((LA16_78 >= '\u0000' && LA16_78 <= 'Z')||LA16_78=='\\'||(LA16_78 >= '^' && LA16_78 <= '\uFFFF')) ) {s = 78;}

                        else if ( (LA16_78=='[') ) {s = 80;}

                        if ( s>=0 ) return s;
                        break;

                    case 4 : 
                        int LA16_0 = input.LA(1);

                        s = -1;
                        if ( (LA16_0=='(') ) {s = 1;}

                        else if ( (LA16_0==')') ) {s = 2;}

                        else if ( (LA16_0==',') ) {s = 3;}

                        else if ( (LA16_0==':') ) {s = 4;}

                        else if ( (LA16_0==';') ) {s = 5;}

                        else if ( (LA16_0=='<') ) {s = 6;}

                        else if ( (LA16_0=='=') ) {s = 7;}

                        else if ( (LA16_0=='>') ) {s = 8;}

                        else if ( (LA16_0=='*') ) {s = 9;}

                        else if ( (LA16_0=='P'||LA16_0=='p') ) {s = 10;}

                        else if ( (LA16_0=='D'||LA16_0=='d') ) {s = 11;}

                        else if ( (LA16_0=='M'||LA16_0=='m') ) {s = 12;}

                        else if ( (LA16_0=='T'||LA16_0=='t') ) {s = 13;}

                        else if ( (LA16_0=='U'||LA16_0=='u') ) {s = 14;}

                        else if ( (LA16_0=='N'||LA16_0=='n') ) {s = 15;}

                        else if ( (LA16_0=='F'||LA16_0=='f') ) {s = 16;}

                        else if ( (LA16_0=='S'||LA16_0=='s') ) {s = 17;}

                        else if ( (LA16_0=='G'||LA16_0=='g') ) {s = 18;}

                        else if ( (LA16_0=='#') ) {s = 19;}

                        else if ( (LA16_0=='B'||LA16_0=='b') ) {s = 20;}

                        else if ( (LA16_0=='E'||LA16_0=='e') ) {s = 21;}

                        else if ( ((LA16_0 >= '-' && LA16_0 <= '9')||LA16_0=='?'||LA16_0=='A'||LA16_0=='C'||(LA16_0 >= 'H' && LA16_0 <= 'L')||LA16_0=='O'||(LA16_0 >= 'Q' && LA16_0 <= 'R')||(LA16_0 >= 'V' && LA16_0 <= 'Z')||LA16_0=='\\'||LA16_0=='_'||LA16_0=='a'||LA16_0=='c'||(LA16_0 >= 'h' && LA16_0 <= 'l')||LA16_0=='o'||(LA16_0 >= 'q' && LA16_0 <= 'r')||(LA16_0 >= 'v' && LA16_0 <= 'z')) ) {s = 22;}

                        else if ( (LA16_0=='\"') ) {s = 23;}

                        else if ( (LA16_0=='{') ) {s = 24;}

                        else if ( (LA16_0=='[') ) {s = 25;}

                        else if ( ((LA16_0 >= '\t' && LA16_0 <= '\n')||LA16_0=='\r'||LA16_0==' ') ) {s = 26;}

                        else if ( ((LA16_0 >= '\u0000' && LA16_0 <= '\b')||(LA16_0 >= '\u000B' && LA16_0 <= '\f')||(LA16_0 >= '\u000E' && LA16_0 <= '\u001F')||LA16_0=='!'||(LA16_0 >= '$' && LA16_0 <= '\'')||LA16_0=='+'||LA16_0=='@'||(LA16_0 >= ']' && LA16_0 <= '^')||LA16_0=='`'||(LA16_0 >= '|' && LA16_0 <= '\uFFFF')) ) {s = 27;}

                        if ( s>=0 ) return s;
                        break;
            }
            NoViableAltException nvae =
                new NoViableAltException(getDescription(), 16, _s, input);
            error(nvae);
            throw nvae;
        }

    }
 

}
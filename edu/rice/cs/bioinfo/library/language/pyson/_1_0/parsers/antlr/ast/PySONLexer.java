// $ANTLR 3.4 D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g 2011-10-12 18:32:28

package edu.rice.cs.bioinfo.library.language.pyson._1_0.parsers.antlr.ast;


import org.antlr.runtime.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;

@SuppressWarnings({"all", "warnings", "unchecked"})
public class PySONLexer extends Lexer {
    public static final int EOF=-1;
    public static final int T__24=24;
    public static final int T__25=25;
    public static final int T__26=26;
    public static final int T__27=27;
    public static final int T__28=28;
    public static final int BEGIN=4;
    public static final int DEFAULT_INDICATOR=5;
    public static final int ELSE=6;
    public static final int END=7;
    public static final int ID=8;
    public static final int ID_SET=9;
    public static final int NESTED_ML_COMMENT=10;
    public static final int NETWORK=11;
    public static final int NETWORKS=12;
    public static final int PHYLONET=13;
    public static final int QUOTE=14;
    public static final int ROOTAGE_QUALIFIER=15;
    public static final int START=16;
    public static final int TAXA_MAP=17;
    public static final int TAXON_SET_LIST=18;
    public static final int TRANSLATE=19;
    public static final int TREE=20;
    public static final int TREES=21;
    public static final int UTREE=22;
    public static final int WS=23;

    @Override   
    public void displayRecognitionError(String[] tokenNames, RecognitionException e) 
    {
    // dont show erors on std error
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
    public String getGrammarFileName() { return "D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g"; }

    // $ANTLR start "T__24"
    public final void mT__24() throws RecognitionException {
        try {
            int _type = T__24;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:13:7: ( '(' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:13:9: '('
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
    // $ANTLR end "T__24"

    // $ANTLR start "T__25"
    public final void mT__25() throws RecognitionException {
        try {
            int _type = T__25;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:14:7: ( ')' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:14:9: ')'
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
    // $ANTLR end "T__25"

    // $ANTLR start "T__26"
    public final void mT__26() throws RecognitionException {
        try {
            int _type = T__26;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:15:7: ( ',' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:15:9: ','
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
    // $ANTLR end "T__26"

    // $ANTLR start "T__27"
    public final void mT__27() throws RecognitionException {
        try {
            int _type = T__27;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:16:7: ( ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:16:9: ';'
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
    // $ANTLR end "T__27"

    // $ANTLR start "T__28"
    public final void mT__28() throws RecognitionException {
        try {
            int _type = T__28;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:17:7: ( '=' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:17:9: '='
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
    // $ANTLR end "T__28"

    // $ANTLR start "DEFAULT_INDICATOR"
    public final void mDEFAULT_INDICATOR() throws RecognitionException {
        try {
            int _type = DEFAULT_INDICATOR;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:125:2: ( '*' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:125:4: '*'
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:127:9: ( ( 'P' | 'p' ) ( 'H' | 'h' ) ( 'Y' | 'y' ) ( 'L' | 'l' ) ( 'O' | 'o' ) ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:127:11: ( 'P' | 'p' ) ( 'H' | 'h' ) ( 'Y' | 'y' ) ( 'L' | 'l' ) ( 'O' | 'o' ) ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' )
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

    // $ANTLR start "TREE"
    public final void mTREE() throws RecognitionException {
        try {
            int _type = TREE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:129:7: ( ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:129:10: ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:131:8: ( ( 'U' | 'u' ) ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:131:11: ( 'U' | 'u' ) ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:133:9: ( ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:133:11: ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' )
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

    // $ANTLR start "START"
    public final void mSTART() throws RecognitionException {
        try {
            int _type = START;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:136:2: ( '#' ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'X' | 'x' ) ( 'U' | 'u' ) ( 'S' | 's' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:136:4: '#' ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'X' | 'x' ) ( 'U' | 'u' ) ( 'S' | 's' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:138:8: ( ( 'B' | 'b' ) ( 'E' | 'e' ) ( 'G' | 'g' ) ( 'I' | 'i' ) ( 'N' | 'n' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:138:10: ( 'B' | 'b' ) ( 'E' | 'e' ) ( 'G' | 'g' ) ( 'I' | 'i' ) ( 'N' | 'n' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:141:2: ( ( 'T' | 't' ) ( 'R' | 'r' ) ( 'A' | 'a' ) ( 'N' | 'n' ) ( 'S' | 's' ) ( 'L' | 'l' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'E' | 'e' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:141:4: ( 'T' | 't' ) ( 'R' | 'r' ) ( 'A' | 'a' ) ( 'N' | 'n' ) ( 'S' | 's' ) ( 'L' | 'l' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'E' | 'e' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:144:2: ( ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' ) ( 'S' | 's' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:144:4: ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' ) ( 'S' | 's' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:146:7: ( ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) ( 'S' | 's' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:146:9: ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) ( 'S' | 's' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:148:5: ( ( 'E' | 'e' ) ( 'N' | 'n' ) ( 'D' | 'd' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:148:7: ( 'E' | 'e' ) ( 'N' | 'n' ) ( 'D' | 'd' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:150:4: ( ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' | '-' )+ )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:150:6: ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' | '-' )+
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:150:6: ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' | '-' )+
            int cnt1=0;
            loop1:
            do {
                int alt1=2;
                int LA1_0 = input.LA(1);

                if ( (LA1_0=='-'||(LA1_0 >= '0' && LA1_0 <= '9')||(LA1_0 >= 'A' && LA1_0 <= 'Z')||LA1_0=='_'||(LA1_0 >= 'a' && LA1_0 <= 'z')) ) {
                    alt1=1;
                }


                switch (alt1) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
            	    {
            	    if ( input.LA(1)=='-'||(input.LA(1) >= '0' && input.LA(1) <= '9')||(input.LA(1) >= 'A' && input.LA(1) <= 'Z')||input.LA(1)=='_'||(input.LA(1) >= 'a' && input.LA(1) <= 'z') ) {
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:152:7: ( '\"' (~ ( '\"' ) )* '\"' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:152:9: '\"' (~ ( '\"' ) )* '\"'
            {
            match('\"'); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:152:13: (~ ( '\"' ) )*
            loop2:
            do {
                int alt2=2;
                int LA2_0 = input.LA(1);

                if ( ((LA2_0 >= '\u0000' && LA2_0 <= '!')||(LA2_0 >= '#' && LA2_0 <= '\uFFFF')) ) {
                    alt2=1;
                }


                switch (alt2) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
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

    // $ANTLR start "TAXA_MAP"
    public final void mTAXA_MAP() throws RecognitionException {
        try {
            int _type = TAXA_MAP;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:9: ( ( ( WS )* ID ( WS )* ':' ( WS )* ID ( ( WS )* ',' ( WS )* ID )* ( WS )* ';' )+ )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:11: ( ( WS )* ID ( WS )* ':' ( WS )* ID ( ( WS )* ',' ( WS )* ID )* ( WS )* ';' )+
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:11: ( ( WS )* ID ( WS )* ':' ( WS )* ID ( ( WS )* ',' ( WS )* ID )* ( WS )* ';' )+
            int cnt10=0;
            loop10:
            do {
                int alt10=2;
                int LA10_0 = input.LA(1);

                if ( ((LA10_0 >= '\t' && LA10_0 <= '\n')||LA10_0=='\r'||LA10_0==' '||LA10_0=='-'||(LA10_0 >= '0' && LA10_0 <= '9')||(LA10_0 >= 'A' && LA10_0 <= 'Z')||LA10_0=='_'||(LA10_0 >= 'a' && LA10_0 <= 'z')) ) {
                    alt10=1;
                }


                switch (alt10) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:12: ( WS )* ID ( WS )* ':' ( WS )* ID ( ( WS )* ',' ( WS )* ID )* ( WS )* ';'
            	    {
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:12: ( WS )*
            	    loop3:
            	    do {
            	        int alt3=2;
            	        int LA3_0 = input.LA(1);

            	        if ( ((LA3_0 >= '\t' && LA3_0 <= '\n')||LA3_0=='\r'||LA3_0==' ') ) {
            	            alt3=1;
            	        }


            	        switch (alt3) {
            	    	case 1 :
            	    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:12: WS
            	    	    {
            	    	    mWS(); 


            	    	    }
            	    	    break;

            	    	default :
            	    	    break loop3;
            	        }
            	    } while (true);


            	    mID(); 


            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:19: ( WS )*
            	    loop4:
            	    do {
            	        int alt4=2;
            	        int LA4_0 = input.LA(1);

            	        if ( ((LA4_0 >= '\t' && LA4_0 <= '\n')||LA4_0=='\r'||LA4_0==' ') ) {
            	            alt4=1;
            	        }


            	        switch (alt4) {
            	    	case 1 :
            	    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:19: WS
            	    	    {
            	    	    mWS(); 


            	    	    }
            	    	    break;

            	    	default :
            	    	    break loop4;
            	        }
            	    } while (true);


            	    match(':'); 

            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:27: ( WS )*
            	    loop5:
            	    do {
            	        int alt5=2;
            	        int LA5_0 = input.LA(1);

            	        if ( ((LA5_0 >= '\t' && LA5_0 <= '\n')||LA5_0=='\r'||LA5_0==' ') ) {
            	            alt5=1;
            	        }


            	        switch (alt5) {
            	    	case 1 :
            	    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:27: WS
            	    	    {
            	    	    mWS(); 


            	    	    }
            	    	    break;

            	    	default :
            	    	    break loop5;
            	        }
            	    } while (true);


            	    mID(); 


            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:34: ( ( WS )* ',' ( WS )* ID )*
            	    loop8:
            	    do {
            	        int alt8=2;
            	        alt8 = dfa8.predict(input);
            	        switch (alt8) {
            	    	case 1 :
            	    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:35: ( WS )* ',' ( WS )* ID
            	    	    {
            	    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:35: ( WS )*
            	    	    loop6:
            	    	    do {
            	    	        int alt6=2;
            	    	        int LA6_0 = input.LA(1);

            	    	        if ( ((LA6_0 >= '\t' && LA6_0 <= '\n')||LA6_0=='\r'||LA6_0==' ') ) {
            	    	            alt6=1;
            	    	        }


            	    	        switch (alt6) {
            	    	    	case 1 :
            	    	    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:35: WS
            	    	    	    {
            	    	    	    mWS(); 


            	    	    	    }
            	    	    	    break;

            	    	    	default :
            	    	    	    break loop6;
            	    	        }
            	    	    } while (true);


            	    	    match(','); 

            	    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:43: ( WS )*
            	    	    loop7:
            	    	    do {
            	    	        int alt7=2;
            	    	        int LA7_0 = input.LA(1);

            	    	        if ( ((LA7_0 >= '\t' && LA7_0 <= '\n')||LA7_0=='\r'||LA7_0==' ') ) {
            	    	            alt7=1;
            	    	        }


            	    	        switch (alt7) {
            	    	    	case 1 :
            	    	    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:43: WS
            	    	    	    {
            	    	    	    mWS(); 


            	    	    	    }
            	    	    	    break;

            	    	    	default :
            	    	    	    break loop7;
            	    	        }
            	    	    } while (true);


            	    	    mID(); 


            	    	    }
            	    	    break;

            	    	default :
            	    	    break loop8;
            	        }
            	    } while (true);


            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:52: ( WS )*
            	    loop9:
            	    do {
            	        int alt9=2;
            	        int LA9_0 = input.LA(1);

            	        if ( ((LA9_0 >= '\t' && LA9_0 <= '\n')||LA9_0=='\r'||LA9_0==' ') ) {
            	            alt9=1;
            	        }


            	        switch (alt9) {
            	    	case 1 :
            	    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:52: WS
            	    	    {
            	    	    mWS(); 


            	    	    }
            	    	    break;

            	    	default :
            	    	    break loop9;
            	        }
            	    } while (true);


            	    match(';'); 

            	    }
            	    break;

            	default :
            	    if ( cnt10 >= 1 ) break loop10;
                        EarlyExitException eee =
                            new EarlyExitException(10, input);
                        throw eee;
                }
                cnt10++;
            } while (true);


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "TAXA_MAP"

    // $ANTLR start "ID_SET"
    public final void mID_SET() throws RecognitionException {
        try {
            int _type = ID_SET;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:157:2: ( '{' ( ( WS )* ID ( WS )* ',' )* ( WS )* ID ( WS )* '}' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:157:4: '{' ( ( WS )* ID ( WS )* ',' )* ( WS )* ID ( WS )* '}'
            {
            match('{'); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:157:8: ( ( WS )* ID ( WS )* ',' )*
            loop13:
            do {
                int alt13=2;
                alt13 = dfa13.predict(input);
                switch (alt13) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:157:9: ( WS )* ID ( WS )* ','
            	    {
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:157:9: ( WS )*
            	    loop11:
            	    do {
            	        int alt11=2;
            	        int LA11_0 = input.LA(1);

            	        if ( ((LA11_0 >= '\t' && LA11_0 <= '\n')||LA11_0=='\r'||LA11_0==' ') ) {
            	            alt11=1;
            	        }


            	        switch (alt11) {
            	    	case 1 :
            	    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:157:9: WS
            	    	    {
            	    	    mWS(); 


            	    	    }
            	    	    break;

            	    	default :
            	    	    break loop11;
            	        }
            	    } while (true);


            	    mID(); 


            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:157:16: ( WS )*
            	    loop12:
            	    do {
            	        int alt12=2;
            	        int LA12_0 = input.LA(1);

            	        if ( ((LA12_0 >= '\t' && LA12_0 <= '\n')||LA12_0=='\r'||LA12_0==' ') ) {
            	            alt12=1;
            	        }


            	        switch (alt12) {
            	    	case 1 :
            	    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:157:16: WS
            	    	    {
            	    	    mWS(); 


            	    	    }
            	    	    break;

            	    	default :
            	    	    break loop12;
            	        }
            	    } while (true);


            	    match(','); 

            	    }
            	    break;

            	default :
            	    break loop13;
                }
            } while (true);


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:158:15: ( WS )*
            loop14:
            do {
                int alt14=2;
                int LA14_0 = input.LA(1);

                if ( ((LA14_0 >= '\t' && LA14_0 <= '\n')||LA14_0=='\r'||LA14_0==' ') ) {
                    alt14=1;
                }


                switch (alt14) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:158:15: WS
            	    {
            	    mWS(); 


            	    }
            	    break;

            	default :
            	    break loop14;
                }
            } while (true);


            mID(); 


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:158:22: ( WS )*
            loop15:
            do {
                int alt15=2;
                int LA15_0 = input.LA(1);

                if ( ((LA15_0 >= '\t' && LA15_0 <= '\n')||LA15_0=='\r'||LA15_0==' ') ) {
                    alt15=1;
                }


                switch (alt15) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:158:22: WS
            	    {
            	    mWS(); 


            	    }
            	    break;

            	default :
            	    break loop15;
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:161:2: ( '(' ( ( WS )* ID_SET ( ( WS )* ',' )? )* ( ( WS )* ID_SET ( WS )* ) ')' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:161:4: '(' ( ( WS )* ID_SET ( ( WS )* ',' )? )* ( ( WS )* ID_SET ( WS )* ) ')'
            {
            match('('); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:161:8: ( ( WS )* ID_SET ( ( WS )* ',' )? )*
            loop19:
            do {
                int alt19=2;
                alt19 = dfa19.predict(input);
                switch (alt19) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:161:9: ( WS )* ID_SET ( ( WS )* ',' )?
            	    {
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:161:9: ( WS )*
            	    loop16:
            	    do {
            	        int alt16=2;
            	        int LA16_0 = input.LA(1);

            	        if ( ((LA16_0 >= '\t' && LA16_0 <= '\n')||LA16_0=='\r'||LA16_0==' ') ) {
            	            alt16=1;
            	        }


            	        switch (alt16) {
            	    	case 1 :
            	    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:161:9: WS
            	    	    {
            	    	    mWS(); 


            	    	    }
            	    	    break;

            	    	default :
            	    	    break loop16;
            	        }
            	    } while (true);


            	    mID_SET(); 


            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:161:20: ( ( WS )* ',' )?
            	    int alt18=2;
            	    alt18 = dfa18.predict(input);
            	    switch (alt18) {
            	        case 1 :
            	            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:161:21: ( WS )* ','
            	            {
            	            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:161:21: ( WS )*
            	            loop17:
            	            do {
            	                int alt17=2;
            	                int LA17_0 = input.LA(1);

            	                if ( ((LA17_0 >= '\t' && LA17_0 <= '\n')||LA17_0=='\r'||LA17_0==' ') ) {
            	                    alt17=1;
            	                }


            	                switch (alt17) {
            	            	case 1 :
            	            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:161:21: WS
            	            	    {
            	            	    mWS(); 


            	            	    }
            	            	    break;

            	            	default :
            	            	    break loop17;
            	                }
            	            } while (true);


            	            match(','); 

            	            }
            	            break;

            	    }


            	    }
            	    break;

            	default :
            	    break loop19;
                }
            } while (true);


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:162:14: ( ( WS )* ID_SET ( WS )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:162:15: ( WS )* ID_SET ( WS )*
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:162:15: ( WS )*
            loop20:
            do {
                int alt20=2;
                int LA20_0 = input.LA(1);

                if ( ((LA20_0 >= '\t' && LA20_0 <= '\n')||LA20_0=='\r'||LA20_0==' ') ) {
                    alt20=1;
                }


                switch (alt20) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:162:15: WS
            	    {
            	    mWS(); 


            	    }
            	    break;

            	default :
            	    break loop20;
                }
            } while (true);


            mID_SET(); 


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:162:26: ( WS )*
            loop21:
            do {
                int alt21=2;
                int LA21_0 = input.LA(1);

                if ( ((LA21_0 >= '\t' && LA21_0 <= '\n')||LA21_0=='\r'||LA21_0==' ') ) {
                    alt21=1;
                }


                switch (alt21) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:162:26: WS
            	    {
            	    mWS(); 


            	    }
            	    break;

            	default :
            	    break loop21;
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

    // $ANTLR start "ROOTAGE_QUALIFIER"
    public final void mROOTAGE_QUALIFIER() throws RecognitionException {
        try {
            int _type = ROOTAGE_QUALIFIER;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:165:2: ( '[' '&' ( 'r' | 'R' | 'u' | 'U' ) ']' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:165:4: '[' '&' ( 'r' | 'R' | 'u' | 'U' ) ']'
            {
            match('['); 

            match('&'); 

            if ( input.LA(1)=='R'||input.LA(1)=='U'||input.LA(1)=='r'||input.LA(1)=='u' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            match(']'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "ROOTAGE_QUALIFIER"

    // $ANTLR start "NESTED_ML_COMMENT"
    public final void mNESTED_ML_COMMENT() throws RecognitionException {
        try {
            int _type = NESTED_ML_COMMENT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:168:3: ( '[' (~ ( '[' | ']' ) )* ']' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:168:6: '[' (~ ( '[' | ']' ) )* ']'
            {
            match('['); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:168:10: (~ ( '[' | ']' ) )*
            loop22:
            do {
                int alt22=2;
                int LA22_0 = input.LA(1);

                if ( ((LA22_0 >= '\u0000' && LA22_0 <= 'Z')||LA22_0=='\\'||(LA22_0 >= '^' && LA22_0 <= '\uFFFF')) ) {
                    alt22=1;
                }


                switch (alt22) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
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
            	    break loop22;
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:171:5: ( ( ' ' | '\\t' | '\\r' | '\\n' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:171:9: ( ' ' | '\\t' | '\\r' | '\\n' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:180:6: ( . )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:180:8: .
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
        // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:8: ( T__24 | T__25 | T__26 | T__27 | T__28 | DEFAULT_INDICATOR | PHYLONET | TREE | UTREE | NETWORK | START | BEGIN | TRANSLATE | NETWORKS | TREES | END | ID | QUOTE | TAXA_MAP | ID_SET | TAXON_SET_LIST | ROOTAGE_QUALIFIER | NESTED_ML_COMMENT | WS | ELSE )
        int alt23=25;
        alt23 = dfa23.predict(input);
        switch (alt23) {
            case 1 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:10: T__24
                {
                mT__24(); 


                }
                break;
            case 2 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:16: T__25
                {
                mT__25(); 


                }
                break;
            case 3 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:22: T__26
                {
                mT__26(); 


                }
                break;
            case 4 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:28: T__27
                {
                mT__27(); 


                }
                break;
            case 5 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:34: T__28
                {
                mT__28(); 


                }
                break;
            case 6 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:40: DEFAULT_INDICATOR
                {
                mDEFAULT_INDICATOR(); 


                }
                break;
            case 7 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:58: PHYLONET
                {
                mPHYLONET(); 


                }
                break;
            case 8 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:67: TREE
                {
                mTREE(); 


                }
                break;
            case 9 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:72: UTREE
                {
                mUTREE(); 


                }
                break;
            case 10 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:78: NETWORK
                {
                mNETWORK(); 


                }
                break;
            case 11 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:86: START
                {
                mSTART(); 


                }
                break;
            case 12 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:92: BEGIN
                {
                mBEGIN(); 


                }
                break;
            case 13 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:98: TRANSLATE
                {
                mTRANSLATE(); 


                }
                break;
            case 14 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:108: NETWORKS
                {
                mNETWORKS(); 


                }
                break;
            case 15 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:117: TREES
                {
                mTREES(); 


                }
                break;
            case 16 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:123: END
                {
                mEND(); 


                }
                break;
            case 17 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:127: ID
                {
                mID(); 


                }
                break;
            case 18 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:130: QUOTE
                {
                mQUOTE(); 


                }
                break;
            case 19 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:136: TAXA_MAP
                {
                mTAXA_MAP(); 


                }
                break;
            case 20 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:145: ID_SET
                {
                mID_SET(); 


                }
                break;
            case 21 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:152: TAXON_SET_LIST
                {
                mTAXON_SET_LIST(); 


                }
                break;
            case 22 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:167: ROOTAGE_QUALIFIER
                {
                mROOTAGE_QUALIFIER(); 


                }
                break;
            case 23 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:185: NESTED_ML_COMMENT
                {
                mNESTED_ML_COMMENT(); 


                }
                break;
            case 24 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:203: WS
                {
                mWS(); 


                }
                break;
            case 25 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:206: ELSE
                {
                mELSE(); 


                }
                break;

        }

    }


    protected DFA8 dfa8 = new DFA8(this);
    protected DFA13 dfa13 = new DFA13(this);
    protected DFA19 dfa19 = new DFA19(this);
    protected DFA18 dfa18 = new DFA18(this);
    protected DFA23 dfa23 = new DFA23(this);
    static final String DFA8_eotS =
        "\4\uffff";
    static final String DFA8_eofS =
        "\4\uffff";
    static final String DFA8_minS =
        "\2\11\2\uffff";
    static final String DFA8_maxS =
        "\2\73\2\uffff";
    static final String DFA8_acceptS =
        "\2\uffff\1\2\1\1";
    static final String DFA8_specialS =
        "\4\uffff}>";
    static final String[] DFA8_transitionS = {
            "\2\1\2\uffff\1\1\22\uffff\1\1\13\uffff\1\3\16\uffff\1\2",
            "\2\1\2\uffff\1\1\22\uffff\1\1\13\uffff\1\3\16\uffff\1\2",
            "",
            ""
    };

    static final short[] DFA8_eot = DFA.unpackEncodedString(DFA8_eotS);
    static final short[] DFA8_eof = DFA.unpackEncodedString(DFA8_eofS);
    static final char[] DFA8_min = DFA.unpackEncodedStringToUnsignedChars(DFA8_minS);
    static final char[] DFA8_max = DFA.unpackEncodedStringToUnsignedChars(DFA8_maxS);
    static final short[] DFA8_accept = DFA.unpackEncodedString(DFA8_acceptS);
    static final short[] DFA8_special = DFA.unpackEncodedString(DFA8_specialS);
    static final short[][] DFA8_transition;

    static {
        int numStates = DFA8_transitionS.length;
        DFA8_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA8_transition[i] = DFA.unpackEncodedString(DFA8_transitionS[i]);
        }
    }

    class DFA8 extends DFA {

        public DFA8(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 8;
            this.eot = DFA8_eot;
            this.eof = DFA8_eof;
            this.min = DFA8_min;
            this.max = DFA8_max;
            this.accept = DFA8_accept;
            this.special = DFA8_special;
            this.transition = DFA8_transition;
        }
        public String getDescription() {
            return "()* loopback of 154:34: ( ( WS )* ',' ( WS )* ID )*";
        }
    }
    static final String DFA13_eotS =
        "\6\uffff";
    static final String DFA13_eofS =
        "\6\uffff";
    static final String DFA13_minS =
        "\4\11\2\uffff";
    static final String DFA13_maxS =
        "\2\172\2\175\2\uffff";
    static final String DFA13_acceptS =
        "\4\uffff\1\2\1\1";
    static final String DFA13_specialS =
        "\6\uffff}>";
    static final String[] DFA13_transitionS = {
            "\2\1\2\uffff\1\1\22\uffff\1\1\14\uffff\1\2\2\uffff\12\2\7\uffff"+
            "\32\2\4\uffff\1\2\1\uffff\32\2",
            "\2\1\2\uffff\1\1\22\uffff\1\1\14\uffff\1\2\2\uffff\12\2\7\uffff"+
            "\32\2\4\uffff\1\2\1\uffff\32\2",
            "\2\3\2\uffff\1\3\22\uffff\1\3\13\uffff\1\5\1\2\2\uffff\12\2"+
            "\7\uffff\32\2\4\uffff\1\2\1\uffff\32\2\2\uffff\1\4",
            "\2\3\2\uffff\1\3\22\uffff\1\3\13\uffff\1\5\120\uffff\1\4",
            "",
            ""
    };

    static final short[] DFA13_eot = DFA.unpackEncodedString(DFA13_eotS);
    static final short[] DFA13_eof = DFA.unpackEncodedString(DFA13_eofS);
    static final char[] DFA13_min = DFA.unpackEncodedStringToUnsignedChars(DFA13_minS);
    static final char[] DFA13_max = DFA.unpackEncodedStringToUnsignedChars(DFA13_maxS);
    static final short[] DFA13_accept = DFA.unpackEncodedString(DFA13_acceptS);
    static final short[] DFA13_special = DFA.unpackEncodedString(DFA13_specialS);
    static final short[][] DFA13_transition;

    static {
        int numStates = DFA13_transitionS.length;
        DFA13_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA13_transition[i] = DFA.unpackEncodedString(DFA13_transitionS[i]);
        }
    }

    class DFA13 extends DFA {

        public DFA13(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 13;
            this.eot = DFA13_eot;
            this.eof = DFA13_eof;
            this.min = DFA13_min;
            this.max = DFA13_max;
            this.accept = DFA13_accept;
            this.special = DFA13_special;
            this.transition = DFA13_transition;
        }
        public String getDescription() {
            return "()* loopback of 157:8: ( ( WS )* ID ( WS )* ',' )*";
        }
    }
    static final String DFA19_eotS =
        "\13\uffff";
    static final String DFA19_eofS =
        "\13\uffff";
    static final String DFA19_minS =
        "\11\11\2\uffff";
    static final String DFA19_maxS =
        "\2\173\2\172\2\175\1\172\2\173\2\uffff";
    static final String DFA19_acceptS =
        "\11\uffff\1\2\1\1";
    static final String DFA19_specialS =
        "\13\uffff}>";
    static final String[] DFA19_transitionS = {
            "\2\1\2\uffff\1\1\22\uffff\1\1\132\uffff\1\2",
            "\2\1\2\uffff\1\1\22\uffff\1\1\132\uffff\1\2",
            "\2\3\2\uffff\1\3\22\uffff\1\3\14\uffff\1\4\2\uffff\12\4\7\uffff"+
            "\32\4\4\uffff\1\4\1\uffff\32\4",
            "\2\3\2\uffff\1\3\22\uffff\1\3\14\uffff\1\4\2\uffff\12\4\7\uffff"+
            "\32\4\4\uffff\1\4\1\uffff\32\4",
            "\2\5\2\uffff\1\5\22\uffff\1\5\13\uffff\1\6\1\4\2\uffff\12\4"+
            "\7\uffff\32\4\4\uffff\1\4\1\uffff\32\4\2\uffff\1\7",
            "\2\5\2\uffff\1\5\22\uffff\1\5\13\uffff\1\6\120\uffff\1\7",
            "\2\3\2\uffff\1\3\22\uffff\1\3\14\uffff\1\4\2\uffff\12\4\7\uffff"+
            "\32\4\4\uffff\1\4\1\uffff\32\4",
            "\2\10\2\uffff\1\10\22\uffff\1\10\10\uffff\1\11\2\uffff\1\12"+
            "\116\uffff\1\12",
            "\2\10\2\uffff\1\10\22\uffff\1\10\10\uffff\1\11\2\uffff\1\12"+
            "\116\uffff\1\12",
            "",
            ""
    };

    static final short[] DFA19_eot = DFA.unpackEncodedString(DFA19_eotS);
    static final short[] DFA19_eof = DFA.unpackEncodedString(DFA19_eofS);
    static final char[] DFA19_min = DFA.unpackEncodedStringToUnsignedChars(DFA19_minS);
    static final char[] DFA19_max = DFA.unpackEncodedStringToUnsignedChars(DFA19_maxS);
    static final short[] DFA19_accept = DFA.unpackEncodedString(DFA19_acceptS);
    static final short[] DFA19_special = DFA.unpackEncodedString(DFA19_specialS);
    static final short[][] DFA19_transition;

    static {
        int numStates = DFA19_transitionS.length;
        DFA19_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA19_transition[i] = DFA.unpackEncodedString(DFA19_transitionS[i]);
        }
    }

    class DFA19 extends DFA {

        public DFA19(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 19;
            this.eot = DFA19_eot;
            this.eof = DFA19_eof;
            this.min = DFA19_min;
            this.max = DFA19_max;
            this.accept = DFA19_accept;
            this.special = DFA19_special;
            this.transition = DFA19_transition;
        }
        public String getDescription() {
            return "()* loopback of 161:8: ( ( WS )* ID_SET ( ( WS )* ',' )? )*";
        }
    }
    static final String DFA18_eotS =
        "\4\uffff";
    static final String DFA18_eofS =
        "\4\uffff";
    static final String DFA18_minS =
        "\2\11\2\uffff";
    static final String DFA18_maxS =
        "\2\173\2\uffff";
    static final String DFA18_acceptS =
        "\2\uffff\1\1\1\2";
    static final String DFA18_specialS =
        "\4\uffff}>";
    static final String[] DFA18_transitionS = {
            "\2\1\2\uffff\1\1\22\uffff\1\1\13\uffff\1\2\116\uffff\1\3",
            "\2\1\2\uffff\1\1\22\uffff\1\1\13\uffff\1\2\116\uffff\1\3",
            "",
            ""
    };

    static final short[] DFA18_eot = DFA.unpackEncodedString(DFA18_eotS);
    static final short[] DFA18_eof = DFA.unpackEncodedString(DFA18_eofS);
    static final char[] DFA18_min = DFA.unpackEncodedStringToUnsignedChars(DFA18_minS);
    static final char[] DFA18_max = DFA.unpackEncodedStringToUnsignedChars(DFA18_maxS);
    static final short[] DFA18_accept = DFA.unpackEncodedString(DFA18_acceptS);
    static final short[] DFA18_special = DFA.unpackEncodedString(DFA18_specialS);
    static final short[][] DFA18_transition;

    static {
        int numStates = DFA18_transitionS.length;
        DFA18_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA18_transition[i] = DFA.unpackEncodedString(DFA18_transitionS[i]);
        }
    }

    class DFA18 extends DFA {

        public DFA18(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 18;
            this.eot = DFA18_eot;
            this.eof = DFA18_eof;
            this.min = DFA18_min;
            this.max = DFA18_max;
            this.accept = DFA18_accept;
            this.special = DFA18_special;
            this.transition = DFA18_transition;
        }
        public String getDescription() {
            return "161:20: ( ( WS )* ',' )?";
        }
    }
    static final String DFA23_eotS =
        "\1\uffff\1\24\5\uffff\4\34\1\23\3\34\1\23\1\46\2\23\10\uffff\1\34"+
        "\1\uffff\1\34\1\uffff\3\34\1\uffff\2\34\5\uffff\6\34\1\70\1\uffff"+
        "\1\34\1\74\4\34\2\uffff\1\34\1\103\1\uffff\1\34\1\105\1\34\1\107"+
        "\1\uffff\1\34\1\uffff\1\34\1\uffff\1\34\1\uffff\2\34\1\116\1\117"+
        "\1\34\1\121\2\uffff\1\122\2\uffff";
    static final String DFA23_eofS =
        "\123\uffff";
    static final String DFA23_minS =
        "\1\0\1\11\5\uffff\4\11\1\116\3\11\1\0\2\11\1\0\10\uffff\1\11\1\uffff"+
        "\1\11\1\uffff\3\11\1\uffff\2\11\3\uffff\1\0\1\uffff\7\11\1\0\6\11"+
        "\2\uffff\2\11\1\uffff\4\11\1\uffff\1\11\1\uffff\1\11\1\uffff\1\11"+
        "\1\uffff\6\11\2\uffff\1\11\2\uffff";
    static final String DFA23_maxS =
        "\1\uffff\1\173\5\uffff\4\172\1\156\3\172\1\uffff\2\172\1\uffff\10"+
        "\uffff\1\172\1\uffff\1\172\1\uffff\3\172\1\uffff\2\172\3\uffff\1"+
        "\uffff\1\uffff\7\172\1\uffff\6\172\2\uffff\2\172\1\uffff\4\172\1"+
        "\uffff\1\172\1\uffff\1\172\1\uffff\1\172\1\uffff\6\172\2\uffff\1"+
        "\172\2\uffff";
    static final String DFA23_acceptS =
        "\2\uffff\1\2\1\3\1\4\1\5\1\6\14\uffff\1\31\1\1\1\25\1\2\1\3\1\4"+
        "\1\5\1\6\1\uffff\1\21\1\uffff\1\23\3\uffff\1\13\2\uffff\1\22\1\30"+
        "\1\24\1\uffff\1\27\16\uffff\1\20\1\26\2\uffff\1\10\4\uffff\1\26"+
        "\1\uffff\1\17\1\uffff\1\11\1\uffff\1\14\6\uffff\1\12\1\7\1\uffff"+
        "\1\16\1\15";
    static final String DFA23_specialS =
        "\1\0\16\uffff\1\3\2\uffff\1\4\25\uffff\1\2\10\uffff\1\1\41\uffff}>";
    static final String[] DFA23_transitionS = {
            "\11\23\2\20\2\23\1\20\22\23\1\20\1\23\1\17\1\13\4\23\1\1\1\2"+
            "\1\6\1\23\1\3\1\16\2\23\12\16\1\23\1\4\1\23\1\5\3\23\1\16\1"+
            "\14\2\16\1\15\10\16\1\12\1\16\1\7\3\16\1\10\1\11\5\16\1\22\3"+
            "\23\1\16\1\23\1\16\1\14\2\16\1\15\10\16\1\12\1\16\1\7\3\16\1"+
            "\10\1\11\5\16\1\21\uff84\23",
            "\2\25\2\uffff\1\25\22\uffff\1\25\132\uffff\1\25",
            "",
            "",
            "",
            "",
            "",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\7\35\1\33\22\35\4\uffff\1\35\1\uffff\7\35\1\33"+
            "\22\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\21\35\1\37\10\35\4\uffff\1\35\1\uffff\21\35\1"+
            "\37\10\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\23\35\1\40\6\35\4\uffff\1\35\1\uffff\23\35\1\40"+
            "\6\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\4\35\1\41\25\35\4\uffff\1\35\1\uffff\4\35\1\41"+
            "\25\35",
            "\1\42\37\uffff\1\42",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\4\35\1\43\25\35\4\uffff\1\35\1\uffff\4\35\1\43"+
            "\25\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\15\35\1\44\14\35\4\uffff\1\35\1\uffff\15\35\1"+
            "\44\14\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\32\35\4\uffff\1\35\1\uffff\32\35",
            "\0\45",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\36\2\uffff\12\36"+
            "\7\uffff\32\36\4\uffff\1\36\1\uffff\32\36",
            "\2\47\2\uffff\1\47\22\uffff\1\47\14\uffff\1\47\2\uffff\12\47"+
            "\7\uffff\32\47\4\uffff\1\47\1\uffff\32\47",
            "\46\51\1\50\64\51\1\uffff\uffa4\51",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\30\35\1\52\1\35\4\uffff\1\35\1\uffff\30\35\1\52"+
            "\1\35",
            "",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\32\35\4\uffff\1\35\1\uffff\32\35",
            "",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\1\54\3\35\1\53\25\35\4\uffff\1\35\1\uffff\1\54"+
            "\3\35\1\53\25\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\21\35\1\55\10\35\4\uffff\1\35\1\uffff\21\35\1"+
            "\55\10\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\23\35\1\56\6\35\4\uffff\1\35\1\uffff\23\35\1\56"+
            "\6\35",
            "",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\6\35\1\57\23\35\4\uffff\1\35\1\uffff\6\35\1\57"+
            "\23\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\3\35\1\60\26\35\4\uffff\1\35\1\uffff\3\35\1\60"+
            "\26\35",
            "",
            "",
            "",
            "\122\51\1\61\2\51\1\61\5\51\1\uffff\26\51\1\61\2\51\1\61\uff8a"+
            "\51",
            "",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\13\35\1\62\16\35\4\uffff\1\35\1\uffff\13\35\1"+
            "\62\16\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\4\35\1\63\25\35\4\uffff\1\35\1\uffff\4\35\1\63"+
            "\25\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\15\35\1\64\14\35\4\uffff\1\35\1\uffff\15\35\1"+
            "\64\14\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\4\35\1\65\25\35\4\uffff\1\35\1\uffff\4\35\1\65"+
            "\25\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\26\35\1\66\3\35\4\uffff\1\35\1\uffff\26\35\1\66"+
            "\3\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\10\35\1\67\21\35\4\uffff\1\35\1\uffff\10\35\1"+
            "\67\21\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\32\35\4\uffff\1\35\1\uffff\32\35",
            "\133\51\1\uffff\1\51\1\71\uffa2\51",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\16\35\1\72\13\35\4\uffff\1\35\1\uffff\16\35\1"+
            "\72\13\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\22\35\1\73\7\35\4\uffff\1\35\1\uffff\22\35\1\73"+
            "\7\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\22\35\1\75\7\35\4\uffff\1\35\1\uffff\22\35\1\75"+
            "\7\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\4\35\1\76\25\35\4\uffff\1\35\1\uffff\4\35\1\76"+
            "\25\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\16\35\1\77\13\35\4\uffff\1\35\1\uffff\16\35\1"+
            "\77\13\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\15\35\1\100\14\35\4\uffff\1\35\1\uffff\15\35\1"+
            "\100\14\35",
            "",
            "",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\15\35\1\102\14\35\4\uffff\1\35\1\uffff\15\35\1"+
            "\102\14\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\32\35\4\uffff\1\35\1\uffff\32\35",
            "",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\13\35\1\104\16\35\4\uffff\1\35\1\uffff\13\35\1"+
            "\104\16\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\32\35\4\uffff\1\35\1\uffff\32\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\21\35\1\106\10\35\4\uffff\1\35\1\uffff\21\35\1"+
            "\106\10\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\32\35\4\uffff\1\35\1\uffff\32\35",
            "",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\4\35\1\110\25\35\4\uffff\1\35\1\uffff\4\35\1\110"+
            "\25\35",
            "",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\1\111\31\35\4\uffff\1\35\1\uffff\1\111\31\35",
            "",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\12\35\1\112\17\35\4\uffff\1\35\1\uffff\12\35\1"+
            "\112\17\35",
            "",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\23\35\1\113\6\35\4\uffff\1\35\1\uffff\23\35\1"+
            "\113\6\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\23\35\1\114\6\35\4\uffff\1\35\1\uffff\23\35\1"+
            "\114\6\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\22\35\1\115\7\35\4\uffff\1\35\1\uffff\22\35\1"+
            "\115\7\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\32\35\4\uffff\1\35\1\uffff\32\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\4\35\1\120\25\35\4\uffff\1\35\1\uffff\4\35\1\120"+
            "\25\35",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\32\35\4\uffff\1\35\1\uffff\32\35",
            "",
            "",
            "\2\36\2\uffff\1\36\22\uffff\1\36\14\uffff\1\35\2\uffff\12\35"+
            "\1\36\6\uffff\32\35\4\uffff\1\35\1\uffff\32\35",
            "",
            ""
    };

    static final short[] DFA23_eot = DFA.unpackEncodedString(DFA23_eotS);
    static final short[] DFA23_eof = DFA.unpackEncodedString(DFA23_eofS);
    static final char[] DFA23_min = DFA.unpackEncodedStringToUnsignedChars(DFA23_minS);
    static final char[] DFA23_max = DFA.unpackEncodedStringToUnsignedChars(DFA23_maxS);
    static final short[] DFA23_accept = DFA.unpackEncodedString(DFA23_acceptS);
    static final short[] DFA23_special = DFA.unpackEncodedString(DFA23_specialS);
    static final short[][] DFA23_transition;

    static {
        int numStates = DFA23_transitionS.length;
        DFA23_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA23_transition[i] = DFA.unpackEncodedString(DFA23_transitionS[i]);
        }
    }

    class DFA23 extends DFA {

        public DFA23(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 23;
            this.eot = DFA23_eot;
            this.eof = DFA23_eof;
            this.min = DFA23_min;
            this.max = DFA23_max;
            this.accept = DFA23_accept;
            this.special = DFA23_special;
            this.transition = DFA23_transition;
        }
        public String getDescription() {
            return "1:1: Tokens : ( T__24 | T__25 | T__26 | T__27 | T__28 | DEFAULT_INDICATOR | PHYLONET | TREE | UTREE | NETWORK | START | BEGIN | TRANSLATE | NETWORKS | TREES | END | ID | QUOTE | TAXA_MAP | ID_SET | TAXON_SET_LIST | ROOTAGE_QUALIFIER | NESTED_ML_COMMENT | WS | ELSE );";
        }
        public int specialStateTransition(int s, IntStream _input) throws NoViableAltException {
            IntStream input = _input;
        	int _s = s;
            switch ( s ) {
                    case 0 : 
                        int LA23_0 = input.LA(1);

                        s = -1;
                        if ( (LA23_0=='(') ) {s = 1;}

                        else if ( (LA23_0==')') ) {s = 2;}

                        else if ( (LA23_0==',') ) {s = 3;}

                        else if ( (LA23_0==';') ) {s = 4;}

                        else if ( (LA23_0=='=') ) {s = 5;}

                        else if ( (LA23_0=='*') ) {s = 6;}

                        else if ( (LA23_0=='P'||LA23_0=='p') ) {s = 7;}

                        else if ( (LA23_0=='T'||LA23_0=='t') ) {s = 8;}

                        else if ( (LA23_0=='U'||LA23_0=='u') ) {s = 9;}

                        else if ( (LA23_0=='N'||LA23_0=='n') ) {s = 10;}

                        else if ( (LA23_0=='#') ) {s = 11;}

                        else if ( (LA23_0=='B'||LA23_0=='b') ) {s = 12;}

                        else if ( (LA23_0=='E'||LA23_0=='e') ) {s = 13;}

                        else if ( (LA23_0=='-'||(LA23_0 >= '0' && LA23_0 <= '9')||LA23_0=='A'||(LA23_0 >= 'C' && LA23_0 <= 'D')||(LA23_0 >= 'F' && LA23_0 <= 'M')||LA23_0=='O'||(LA23_0 >= 'Q' && LA23_0 <= 'S')||(LA23_0 >= 'V' && LA23_0 <= 'Z')||LA23_0=='_'||LA23_0=='a'||(LA23_0 >= 'c' && LA23_0 <= 'd')||(LA23_0 >= 'f' && LA23_0 <= 'm')||LA23_0=='o'||(LA23_0 >= 'q' && LA23_0 <= 's')||(LA23_0 >= 'v' && LA23_0 <= 'z')) ) {s = 14;}

                        else if ( (LA23_0=='\"') ) {s = 15;}

                        else if ( ((LA23_0 >= '\t' && LA23_0 <= '\n')||LA23_0=='\r'||LA23_0==' ') ) {s = 16;}

                        else if ( (LA23_0=='{') ) {s = 17;}

                        else if ( (LA23_0=='[') ) {s = 18;}

                        else if ( ((LA23_0 >= '\u0000' && LA23_0 <= '\b')||(LA23_0 >= '\u000B' && LA23_0 <= '\f')||(LA23_0 >= '\u000E' && LA23_0 <= '\u001F')||LA23_0=='!'||(LA23_0 >= '$' && LA23_0 <= '\'')||LA23_0=='+'||(LA23_0 >= '.' && LA23_0 <= '/')||LA23_0==':'||LA23_0=='<'||(LA23_0 >= '>' && LA23_0 <= '@')||(LA23_0 >= '\\' && LA23_0 <= '^')||LA23_0=='`'||(LA23_0 >= '|' && LA23_0 <= '\uFFFF')) ) {s = 19;}

                        if ( s>=0 ) return s;
                        break;

                    case 1 : 
                        int LA23_49 = input.LA(1);

                        s = -1;
                        if ( (LA23_49==']') ) {s = 57;}

                        else if ( ((LA23_49 >= '\u0000' && LA23_49 <= 'Z')||LA23_49=='\\'||(LA23_49 >= '^' && LA23_49 <= '\uFFFF')) ) {s = 41;}

                        if ( s>=0 ) return s;
                        break;

                    case 2 : 
                        int LA23_40 = input.LA(1);

                        s = -1;
                        if ( (LA23_40=='R'||LA23_40=='U'||LA23_40=='r'||LA23_40=='u') ) {s = 49;}

                        else if ( ((LA23_40 >= '\u0000' && LA23_40 <= 'Q')||(LA23_40 >= 'S' && LA23_40 <= 'T')||(LA23_40 >= 'V' && LA23_40 <= 'Z')||(LA23_40 >= '\\' && LA23_40 <= 'q')||(LA23_40 >= 's' && LA23_40 <= 't')||(LA23_40 >= 'v' && LA23_40 <= '\uFFFF')) ) {s = 41;}

                        if ( s>=0 ) return s;
                        break;

                    case 3 : 
                        int LA23_15 = input.LA(1);

                        s = -1;
                        if ( ((LA23_15 >= '\u0000' && LA23_15 <= '\uFFFF')) ) {s = 37;}

                        else s = 19;

                        if ( s>=0 ) return s;
                        break;

                    case 4 : 
                        int LA23_18 = input.LA(1);

                        s = -1;
                        if ( (LA23_18=='&') ) {s = 40;}

                        else if ( ((LA23_18 >= '\u0000' && LA23_18 <= '%')||(LA23_18 >= '\'' && LA23_18 <= 'Z')||(LA23_18 >= '\\' && LA23_18 <= '\uFFFF')) ) {s = 41;}

                        else s = 19;

                        if ( s>=0 ) return s;
                        break;
            }
            NoViableAltException nvae =
                new NoViableAltException(getDescription(), 23, _s, input);
            error(nvae);
            throw nvae;
        }

    }
 

}
// $ANTLR 3.4 D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g 2011-08-31 18:58:15

package edu.rice.bioinfo.library.language.pyson._1_0.parsers.antlr.ast;


import org.antlr.runtime.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;

@SuppressWarnings({"all", "warnings", "unchecked"})
public class PySONLexer extends Lexer {
    public static final int EOF=-1;
    public static final int T__17=17;
    public static final int T__18=18;
    public static final int T__19=19;
    public static final int BEGIN=4;
    public static final int DEFAULT_INDICATOR=5;
    public static final int END=6;
    public static final int ID=7;
    public static final int NESTED_ML_COMMENT=8;
    public static final int NETWORK=9;
    public static final int NETWORKS=10;
    public static final int START=11;
    public static final int TRANSLATE=12;
    public static final int TREE=13;
    public static final int TREES=14;
    public static final int UTREE=15;
    public static final int WS=16;

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

    // $ANTLR start "T__17"
    public final void mT__17() throws RecognitionException {
        try {
            int _type = T__17;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:6:7: ( ',' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:6:9: ','
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
    // $ANTLR end "T__17"

    // $ANTLR start "T__18"
    public final void mT__18() throws RecognitionException {
        try {
            int _type = T__18;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:7:7: ( ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:7:9: ';'
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
    // $ANTLR end "T__18"

    // $ANTLR start "T__19"
    public final void mT__19() throws RecognitionException {
        try {
            int _type = T__19;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:8:7: ( '=' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:8:9: '='
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
    // $ANTLR end "T__19"

    // $ANTLR start "DEFAULT_INDICATOR"
    public final void mDEFAULT_INDICATOR() throws RecognitionException {
        try {
            int _type = DEFAULT_INDICATOR;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:2: ( '*' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:4: '*'
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

    // $ANTLR start "TREE"
    public final void mTREE() throws RecognitionException {
        try {
            int _type = TREE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:91:7: ( ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:91:10: ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:93:8: ( ( 'U' | 'u' ) ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:93:11: ( 'U' | 'u' ) ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:95:9: ( ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:95:11: ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:98:2: ( '#' ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'X' | 'x' ) ( 'U' | 'u' ) ( 'S' | 's' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:98:4: '#' ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'X' | 'x' ) ( 'U' | 'u' ) ( 'S' | 's' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:100:8: ( ( 'B' | 'b' ) ( 'E' | 'e' ) ( 'G' | 'g' ) ( 'I' | 'i' ) ( 'N' | 'n' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:100:10: ( 'B' | 'b' ) ( 'E' | 'e' ) ( 'G' | 'g' ) ( 'I' | 'i' ) ( 'N' | 'n' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:103:2: ( ( 'T' | 't' ) ( 'R' | 'r' ) ( 'A' | 'a' ) ( 'N' | 'n' ) ( 'S' | 's' ) ( 'L' | 'l' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'E' | 'e' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:103:4: ( 'T' | 't' ) ( 'R' | 'r' ) ( 'A' | 'a' ) ( 'N' | 'n' ) ( 'S' | 's' ) ( 'L' | 'l' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'E' | 'e' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:106:2: ( ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' ) ( 'S' | 's' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:106:4: ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' ) ( 'S' | 's' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:108:7: ( ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) ( 'S' | 's' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:108:9: ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) ( 'S' | 's' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:110:5: ( ( 'E' | 'e' ) ( 'N' | 'n' ) ( 'D' | 'd' ) ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:110:7: ( 'E' | 'e' ) ( 'N' | 'n' ) ( 'D' | 'd' ) ';'
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


            match(';'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "END"

    // $ANTLR start "NESTED_ML_COMMENT"
    public final void mNESTED_ML_COMMENT() throws RecognitionException {
        try {
            int _type = NESTED_ML_COMMENT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:113:3: ( '[' (~ ( '[' | ']' ) | NESTED_ML_COMMENT )* ']' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:113:6: '[' (~ ( '[' | ']' ) | NESTED_ML_COMMENT )* ']'
            {
            match('['); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:113:10: (~ ( '[' | ']' ) | NESTED_ML_COMMENT )*
            loop1:
            do {
                int alt1=3;
                int LA1_0 = input.LA(1);

                if ( ((LA1_0 >= '\u0000' && LA1_0 <= 'Z')||LA1_0=='\\'||(LA1_0 >= '^' && LA1_0 <= '\uFFFF')) ) {
                    alt1=1;
                }
                else if ( (LA1_0=='[') ) {
                    alt1=2;
                }


                switch (alt1) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:113:11: ~ ( '[' | ']' )
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
            	case 2 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:113:26: NESTED_ML_COMMENT
            	    {
            	    mNESTED_ML_COMMENT(); 


            	    }
            	    break;

            	default :
            	    break loop1;
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:115:5: ( ( ' ' | '\\t' | '\\r' | '\\n' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:115:9: ( ' ' | '\\t' | '\\r' | '\\n' )
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

    // $ANTLR start "ID"
    public final void mID() throws RecognitionException {
        try {
            int _type = ID;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:123:5: ( ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) )+ )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:123:7: ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) )+
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:123:7: ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) )+
            int cnt2=0;
            loop2:
            do {
                int alt2=2;
                int LA2_0 = input.LA(1);

                if ( ((LA2_0 >= '0' && LA2_0 <= '9')||(LA2_0 >= 'A' && LA2_0 <= 'Z')||(LA2_0 >= 'a' && LA2_0 <= 'z')) ) {
                    alt2=1;
                }


                switch (alt2) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
            	    {
            	    if ( (input.LA(1) >= '0' && input.LA(1) <= '9')||(input.LA(1) >= 'A' && input.LA(1) <= 'Z')||(input.LA(1) >= 'a' && input.LA(1) <= 'z') ) {
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
            	    if ( cnt2 >= 1 ) break loop2;
                        EarlyExitException eee =
                            new EarlyExitException(2, input);
                        throw eee;
                }
                cnt2++;
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

    public void mTokens() throws RecognitionException {
        // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:8: ( T__17 | T__18 | T__19 | DEFAULT_INDICATOR | TREE | UTREE | NETWORK | START | BEGIN | TRANSLATE | NETWORKS | TREES | END | NESTED_ML_COMMENT | WS | ID )
        int alt3=16;
        switch ( input.LA(1) ) {
        case ',':
            {
            alt3=1;
            }
            break;
        case ';':
            {
            alt3=2;
            }
            break;
        case '=':
            {
            alt3=3;
            }
            break;
        case '*':
            {
            alt3=4;
            }
            break;
        case 'T':
        case 't':
            {
            int LA3_5 = input.LA(2);

            if ( (LA3_5=='R'||LA3_5=='r') ) {
                switch ( input.LA(3) ) {
                case 'E':
                case 'e':
                    {
                    int LA3_19 = input.LA(4);

                    if ( (LA3_19=='E'||LA3_19=='e') ) {
                        switch ( input.LA(5) ) {
                        case 'S':
                        case 's':
                            {
                            int LA3_31 = input.LA(6);

                            if ( ((LA3_31 >= '0' && LA3_31 <= '9')||(LA3_31 >= 'A' && LA3_31 <= 'Z')||(LA3_31 >= 'a' && LA3_31 <= 'z')) ) {
                                alt3=16;
                            }
                            else {
                                alt3=12;
                            }
                            }
                            break;
                        case '0':
                        case '1':
                        case '2':
                        case '3':
                        case '4':
                        case '5':
                        case '6':
                        case '7':
                        case '8':
                        case '9':
                        case 'A':
                        case 'B':
                        case 'C':
                        case 'D':
                        case 'E':
                        case 'F':
                        case 'G':
                        case 'H':
                        case 'I':
                        case 'J':
                        case 'K':
                        case 'L':
                        case 'M':
                        case 'N':
                        case 'O':
                        case 'P':
                        case 'Q':
                        case 'R':
                        case 'T':
                        case 'U':
                        case 'V':
                        case 'W':
                        case 'X':
                        case 'Y':
                        case 'Z':
                        case 'a':
                        case 'b':
                        case 'c':
                        case 'd':
                        case 'e':
                        case 'f':
                        case 'g':
                        case 'h':
                        case 'i':
                        case 'j':
                        case 'k':
                        case 'l':
                        case 'm':
                        case 'n':
                        case 'o':
                        case 'p':
                        case 'q':
                        case 'r':
                        case 't':
                        case 'u':
                        case 'v':
                        case 'w':
                        case 'x':
                        case 'y':
                        case 'z':
                            {
                            alt3=16;
                            }
                            break;
                        default:
                            alt3=5;
                        }

                    }
                    else {
                        alt3=16;
                    }
                    }
                    break;
                case 'A':
                case 'a':
                    {
                    int LA3_20 = input.LA(4);

                    if ( (LA3_20=='N'||LA3_20=='n') ) {
                        int LA3_26 = input.LA(5);

                        if ( (LA3_26=='S'||LA3_26=='s') ) {
                            int LA3_33 = input.LA(6);

                            if ( (LA3_33=='L'||LA3_33=='l') ) {
                                int LA3_38 = input.LA(7);

                                if ( (LA3_38=='A'||LA3_38=='a') ) {
                                    int LA3_42 = input.LA(8);

                                    if ( (LA3_42=='T'||LA3_42=='t') ) {
                                        int LA3_44 = input.LA(9);

                                        if ( (LA3_44=='E'||LA3_44=='e') ) {
                                            int LA3_47 = input.LA(10);

                                            if ( ((LA3_47 >= '0' && LA3_47 <= '9')||(LA3_47 >= 'A' && LA3_47 <= 'Z')||(LA3_47 >= 'a' && LA3_47 <= 'z')) ) {
                                                alt3=16;
                                            }
                                            else {
                                                alt3=10;
                                            }
                                        }
                                        else {
                                            alt3=16;
                                        }
                                    }
                                    else {
                                        alt3=16;
                                    }
                                }
                                else {
                                    alt3=16;
                                }
                            }
                            else {
                                alt3=16;
                            }
                        }
                        else {
                            alt3=16;
                        }
                    }
                    else {
                        alt3=16;
                    }
                    }
                    break;
                default:
                    alt3=16;
                }

            }
            else {
                alt3=16;
            }
            }
            break;
        case 'U':
        case 'u':
            {
            int LA3_6 = input.LA(2);

            if ( (LA3_6=='T'||LA3_6=='t') ) {
                int LA3_15 = input.LA(3);

                if ( (LA3_15=='R'||LA3_15=='r') ) {
                    int LA3_21 = input.LA(4);

                    if ( (LA3_21=='E'||LA3_21=='e') ) {
                        int LA3_27 = input.LA(5);

                        if ( (LA3_27=='E'||LA3_27=='e') ) {
                            int LA3_34 = input.LA(6);

                            if ( ((LA3_34 >= '0' && LA3_34 <= '9')||(LA3_34 >= 'A' && LA3_34 <= 'Z')||(LA3_34 >= 'a' && LA3_34 <= 'z')) ) {
                                alt3=16;
                            }
                            else {
                                alt3=6;
                            }
                        }
                        else {
                            alt3=16;
                        }
                    }
                    else {
                        alt3=16;
                    }
                }
                else {
                    alt3=16;
                }
            }
            else {
                alt3=16;
            }
            }
            break;
        case 'N':
        case 'n':
            {
            int LA3_7 = input.LA(2);

            if ( (LA3_7=='E'||LA3_7=='e') ) {
                int LA3_16 = input.LA(3);

                if ( (LA3_16=='T'||LA3_16=='t') ) {
                    int LA3_22 = input.LA(4);

                    if ( (LA3_22=='W'||LA3_22=='w') ) {
                        int LA3_28 = input.LA(5);

                        if ( (LA3_28=='O'||LA3_28=='o') ) {
                            int LA3_35 = input.LA(6);

                            if ( (LA3_35=='R'||LA3_35=='r') ) {
                                int LA3_40 = input.LA(7);

                                if ( (LA3_40=='K'||LA3_40=='k') ) {
                                    switch ( input.LA(8) ) {
                                    case 'S':
                                    case 's':
                                        {
                                        int LA3_45 = input.LA(9);

                                        if ( ((LA3_45 >= '0' && LA3_45 <= '9')||(LA3_45 >= 'A' && LA3_45 <= 'Z')||(LA3_45 >= 'a' && LA3_45 <= 'z')) ) {
                                            alt3=16;
                                        }
                                        else {
                                            alt3=11;
                                        }
                                        }
                                        break;
                                    case '0':
                                    case '1':
                                    case '2':
                                    case '3':
                                    case '4':
                                    case '5':
                                    case '6':
                                    case '7':
                                    case '8':
                                    case '9':
                                    case 'A':
                                    case 'B':
                                    case 'C':
                                    case 'D':
                                    case 'E':
                                    case 'F':
                                    case 'G':
                                    case 'H':
                                    case 'I':
                                    case 'J':
                                    case 'K':
                                    case 'L':
                                    case 'M':
                                    case 'N':
                                    case 'O':
                                    case 'P':
                                    case 'Q':
                                    case 'R':
                                    case 'T':
                                    case 'U':
                                    case 'V':
                                    case 'W':
                                    case 'X':
                                    case 'Y':
                                    case 'Z':
                                    case 'a':
                                    case 'b':
                                    case 'c':
                                    case 'd':
                                    case 'e':
                                    case 'f':
                                    case 'g':
                                    case 'h':
                                    case 'i':
                                    case 'j':
                                    case 'k':
                                    case 'l':
                                    case 'm':
                                    case 'n':
                                    case 'o':
                                    case 'p':
                                    case 'q':
                                    case 'r':
                                    case 't':
                                    case 'u':
                                    case 'v':
                                    case 'w':
                                    case 'x':
                                    case 'y':
                                    case 'z':
                                        {
                                        alt3=16;
                                        }
                                        break;
                                    default:
                                        alt3=7;
                                    }

                                }
                                else {
                                    alt3=16;
                                }
                            }
                            else {
                                alt3=16;
                            }
                        }
                        else {
                            alt3=16;
                        }
                    }
                    else {
                        alt3=16;
                    }
                }
                else {
                    alt3=16;
                }
            }
            else {
                alt3=16;
            }
            }
            break;
        case '#':
            {
            alt3=8;
            }
            break;
        case 'B':
        case 'b':
            {
            int LA3_9 = input.LA(2);

            if ( (LA3_9=='E'||LA3_9=='e') ) {
                int LA3_17 = input.LA(3);

                if ( (LA3_17=='G'||LA3_17=='g') ) {
                    int LA3_23 = input.LA(4);

                    if ( (LA3_23=='I'||LA3_23=='i') ) {
                        int LA3_29 = input.LA(5);

                        if ( (LA3_29=='N'||LA3_29=='n') ) {
                            int LA3_36 = input.LA(6);

                            if ( ((LA3_36 >= '0' && LA3_36 <= '9')||(LA3_36 >= 'A' && LA3_36 <= 'Z')||(LA3_36 >= 'a' && LA3_36 <= 'z')) ) {
                                alt3=16;
                            }
                            else {
                                alt3=9;
                            }
                        }
                        else {
                            alt3=16;
                        }
                    }
                    else {
                        alt3=16;
                    }
                }
                else {
                    alt3=16;
                }
            }
            else {
                alt3=16;
            }
            }
            break;
        case 'E':
        case 'e':
            {
            int LA3_10 = input.LA(2);

            if ( (LA3_10=='N'||LA3_10=='n') ) {
                int LA3_18 = input.LA(3);

                if ( (LA3_18=='D'||LA3_18=='d') ) {
                    int LA3_24 = input.LA(4);

                    if ( (LA3_24==';') ) {
                        alt3=13;
                    }
                    else {
                        alt3=16;
                    }
                }
                else {
                    alt3=16;
                }
            }
            else {
                alt3=16;
            }
            }
            break;
        case '[':
            {
            alt3=14;
            }
            break;
        case '\t':
        case '\n':
        case '\r':
        case ' ':
            {
            alt3=15;
            }
            break;
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
        case 'A':
        case 'C':
        case 'D':
        case 'F':
        case 'G':
        case 'H':
        case 'I':
        case 'J':
        case 'K':
        case 'L':
        case 'M':
        case 'O':
        case 'P':
        case 'Q':
        case 'R':
        case 'S':
        case 'V':
        case 'W':
        case 'X':
        case 'Y':
        case 'Z':
        case 'a':
        case 'c':
        case 'd':
        case 'f':
        case 'g':
        case 'h':
        case 'i':
        case 'j':
        case 'k':
        case 'l':
        case 'm':
        case 'o':
        case 'p':
        case 'q':
        case 'r':
        case 's':
        case 'v':
        case 'w':
        case 'x':
        case 'y':
        case 'z':
            {
            alt3=16;
            }
            break;
        default:
            NoViableAltException nvae =
                new NoViableAltException("", 3, 0, input);

            throw nvae;

        }

        switch (alt3) {
            case 1 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:10: T__17
                {
                mT__17(); 


                }
                break;
            case 2 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:16: T__18
                {
                mT__18(); 


                }
                break;
            case 3 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:22: T__19
                {
                mT__19(); 


                }
                break;
            case 4 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:28: DEFAULT_INDICATOR
                {
                mDEFAULT_INDICATOR(); 


                }
                break;
            case 5 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:46: TREE
                {
                mTREE(); 


                }
                break;
            case 6 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:51: UTREE
                {
                mUTREE(); 


                }
                break;
            case 7 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:57: NETWORK
                {
                mNETWORK(); 


                }
                break;
            case 8 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:65: START
                {
                mSTART(); 


                }
                break;
            case 9 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:71: BEGIN
                {
                mBEGIN(); 


                }
                break;
            case 10 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:77: TRANSLATE
                {
                mTRANSLATE(); 


                }
                break;
            case 11 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:87: NETWORKS
                {
                mNETWORKS(); 


                }
                break;
            case 12 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:96: TREES
                {
                mTREES(); 


                }
                break;
            case 13 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:102: END
                {
                mEND(); 


                }
                break;
            case 14 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:106: NESTED_ML_COMMENT
                {
                mNESTED_ML_COMMENT(); 


                }
                break;
            case 15 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:124: WS
                {
                mWS(); 


                }
                break;
            case 16 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:127: ID
                {
                mID(); 


                }
                break;

        }

    }


 

}
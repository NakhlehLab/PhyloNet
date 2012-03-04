/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

// $ANTLR 3.4 D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g 2011-08-23 10:27:56

package edu.rice.cs.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast;


import org.antlr.runtime.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;

@SuppressWarnings({"all", "warnings", "unchecked"})
public class ExtendedNewickLexer extends Lexer {
    public static final int EOF=-1;
    public static final int T__10=10;
    public static final int T__11=11;
    public static final int T__12=12;
    public static final int T__13=13;
    public static final int T__14=14;
    public static final int T__15=15;
    public static final int DECIMAL_NUMBER=4;
    public static final int NESTED_ML_COMMENT=5;
    public static final int QUOTED_TEXT=6;
    public static final int ROOTAGE_QUALIFIER=7;
    public static final int UNQUOTED_ALPHA_TEXT=8;
    public static final int WS=9;

    // delegates
    // delegators
    public Lexer[] getDelegates() {
        return new Lexer[] {};
    }

    public ExtendedNewickLexer() {} 
    public ExtendedNewickLexer(CharStream input) {
        this(input, new RecognizerSharedState());
    }
    public ExtendedNewickLexer(CharStream input, RecognizerSharedState state) {
        super(input,state);
    }
    public String getGrammarFileName() { return "D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g"; }

    // $ANTLR start "T__10"
    public final void mT__10() throws RecognitionException {
        try {
            int _type = T__10;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:6:7: ( '#' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:6:9: '#'
            {
            match('#'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__10"

    // $ANTLR start "T__11"
    public final void mT__11() throws RecognitionException {
        try {
            int _type = T__11;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:7:7: ( '(' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:7:9: '('
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
    // $ANTLR end "T__11"

    // $ANTLR start "T__12"
    public final void mT__12() throws RecognitionException {
        try {
            int _type = T__12;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:8:7: ( ')' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:8:9: ')'
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
    // $ANTLR end "T__12"

    // $ANTLR start "T__13"
    public final void mT__13() throws RecognitionException {
        try {
            int _type = T__13;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:9:7: ( ',' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:9:9: ','
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
    // $ANTLR end "T__13"

    // $ANTLR start "T__14"
    public final void mT__14() throws RecognitionException {
        try {
            int _type = T__14;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:10:7: ( ':' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:10:9: ':'
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
    // $ANTLR end "T__14"

    // $ANTLR start "T__15"
    public final void mT__15() throws RecognitionException {
        try {
            int _type = T__15;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:11:7: ( ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:11:9: ';'
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
    // $ANTLR end "T__15"

    // $ANTLR start "ROOTAGE_QUALIFIER"
    public final void mROOTAGE_QUALIFIER() throws RecognitionException {
        try {
            int _type = ROOTAGE_QUALIFIER;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:108:2: ( '[' '&' ( 'r' | 'R' | 'u' | 'U' ) ']' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:108:4: '[' '&' ( 'r' | 'R' | 'u' | 'U' ) ']'
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

    // $ANTLR start "DECIMAL_NUMBER"
    public final void mDECIMAL_NUMBER() throws RecognitionException {
        try {
            int _type = DECIMAL_NUMBER;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:112:16: ( ( ( '0' .. '9' ) | '.' ) ( ( '0' .. '9' ) | ' ' | '.' )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:112:18: ( ( '0' .. '9' ) | '.' ) ( ( '0' .. '9' ) | ' ' | '.' )*
            {
            if ( input.LA(1)=='.'||(input.LA(1) >= '0' && input.LA(1) <= '9') ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:112:37: ( ( '0' .. '9' ) | ' ' | '.' )*
            loop1:
            do {
                int alt1=2;
                int LA1_0 = input.LA(1);

                if ( (LA1_0==' '||LA1_0=='.'||(LA1_0 >= '0' && LA1_0 <= '9')) ) {
                    alt1=1;
                }


                switch (alt1) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:
            	    {
            	    if ( input.LA(1)==' '||input.LA(1)=='.'||(input.LA(1) >= '0' && input.LA(1) <= '9') ) {
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
            	    break loop1;
                }
            } while (true);


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "DECIMAL_NUMBER"

    // $ANTLR start "UNQUOTED_ALPHA_TEXT"
    public final void mUNQUOTED_ALPHA_TEXT() throws RecognitionException {
        try {
            int _type = UNQUOTED_ALPHA_TEXT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:115:2: (~ ( '(' | ')' | '[' | ']' | ':' | ';' | '#' | '\\'' | ',' | ' ' | '.' | ( '0' .. '9' ) | '\\t' | '\\r' | '\\n' ) (~ ( '(' | ')' | '[' | ']' | ':' | ';' | '#' | '\\'' | ',' | '.' | ( '0' .. '9' ) ) )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:115:4: ~ ( '(' | ')' | '[' | ']' | ':' | ';' | '#' | '\\'' | ',' | ' ' | '.' | ( '0' .. '9' ) | '\\t' | '\\r' | '\\n' ) (~ ( '(' | ')' | '[' | ']' | ':' | ';' | '#' | '\\'' | ',' | '.' | ( '0' .. '9' ) ) )*
            {
            if ( (input.LA(1) >= '\u0000' && input.LA(1) <= '\b')||(input.LA(1) >= '\u000B' && input.LA(1) <= '\f')||(input.LA(1) >= '\u000E' && input.LA(1) <= '\u001F')||(input.LA(1) >= '!' && input.LA(1) <= '\"')||(input.LA(1) >= '$' && input.LA(1) <= '&')||(input.LA(1) >= '*' && input.LA(1) <= '+')||input.LA(1)=='-'||input.LA(1)=='/'||(input.LA(1) >= '<' && input.LA(1) <= 'Z')||input.LA(1)=='\\'||(input.LA(1) >= '^' && input.LA(1) <= '\uFFFF') ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:116:32: (~ ( '(' | ')' | '[' | ']' | ':' | ';' | '#' | '\\'' | ',' | '.' | ( '0' .. '9' ) ) )*
            loop2:
            do {
                int alt2=2;
                int LA2_0 = input.LA(1);

                if ( ((LA2_0 >= '\u0000' && LA2_0 <= '\"')||(LA2_0 >= '$' && LA2_0 <= '&')||(LA2_0 >= '*' && LA2_0 <= '+')||LA2_0=='-'||LA2_0=='/'||(LA2_0 >= '<' && LA2_0 <= 'Z')||LA2_0=='\\'||(LA2_0 >= '^' && LA2_0 <= '\uFFFF')) ) {
                    alt2=1;
                }


                switch (alt2) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:
            	    {
            	    if ( (input.LA(1) >= '\u0000' && input.LA(1) <= '\"')||(input.LA(1) >= '$' && input.LA(1) <= '&')||(input.LA(1) >= '*' && input.LA(1) <= '+')||input.LA(1)=='-'||input.LA(1)=='/'||(input.LA(1) >= '<' && input.LA(1) <= 'Z')||input.LA(1)=='\\'||(input.LA(1) >= '^' && input.LA(1) <= '\uFFFF') ) {
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


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "UNQUOTED_ALPHA_TEXT"

    // $ANTLR start "QUOTED_TEXT"
    public final void mQUOTED_TEXT() throws RecognitionException {
        try {
            int _type = QUOTED_TEXT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:119:2: ( '\\'' (~ ( '\\n' | '\\r' | '\\'' ) | ( '\\'\\'' ) )* '\\'' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:119:4: '\\'' (~ ( '\\n' | '\\r' | '\\'' ) | ( '\\'\\'' ) )* '\\''
            {
            match('\''); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:119:9: (~ ( '\\n' | '\\r' | '\\'' ) | ( '\\'\\'' ) )*
            loop3:
            do {
                int alt3=3;
                int LA3_0 = input.LA(1);

                if ( (LA3_0=='\'') ) {
                    int LA3_1 = input.LA(2);

                    if ( (LA3_1=='\'') ) {
                        alt3=2;
                    }


                }
                else if ( ((LA3_0 >= '\u0000' && LA3_0 <= '\t')||(LA3_0 >= '\u000B' && LA3_0 <= '\f')||(LA3_0 >= '\u000E' && LA3_0 <= '&')||(LA3_0 >= '(' && LA3_0 <= '\uFFFF')) ) {
                    alt3=1;
                }


                switch (alt3) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:119:10: ~ ( '\\n' | '\\r' | '\\'' )
            	    {
            	    if ( (input.LA(1) >= '\u0000' && input.LA(1) <= '\t')||(input.LA(1) >= '\u000B' && input.LA(1) <= '\f')||(input.LA(1) >= '\u000E' && input.LA(1) <= '&')||(input.LA(1) >= '(' && input.LA(1) <= '\uFFFF') ) {
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
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:119:34: ( '\\'\\'' )
            	    {
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:119:34: ( '\\'\\'' )
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:119:35: '\\'\\''
            	    {
            	    match("''"); 



            	    }


            	    }
            	    break;

            	default :
            	    break loop3;
                }
            } while (true);


            match('\''); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "QUOTED_TEXT"

    // $ANTLR start "NESTED_ML_COMMENT"
    public final void mNESTED_ML_COMMENT() throws RecognitionException {
        try {
            int _type = NESTED_ML_COMMENT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:122:3: ( '[' (~ ( '[' | ']' ) | NESTED_ML_COMMENT )* ']' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:122:6: '[' (~ ( '[' | ']' ) | NESTED_ML_COMMENT )* ']'
            {
            match('['); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:122:10: (~ ( '[' | ']' ) | NESTED_ML_COMMENT )*
            loop4:
            do {
                int alt4=3;
                int LA4_0 = input.LA(1);

                if ( ((LA4_0 >= '\u0000' && LA4_0 <= 'Z')||LA4_0=='\\'||(LA4_0 >= '^' && LA4_0 <= '\uFFFF')) ) {
                    alt4=1;
                }
                else if ( (LA4_0=='[') ) {
                    alt4=2;
                }


                switch (alt4) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:122:11: ~ ( '[' | ']' )
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
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:122:26: NESTED_ML_COMMENT
            	    {
            	    mNESTED_ML_COMMENT(); 


            	    }
            	    break;

            	default :
            	    break loop4;
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:124:5: ( ( ' ' | '\\t' | '\\r' | '\\n' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:124:9: ( ' ' | '\\t' | '\\r' | '\\n' )
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

    public void mTokens() throws RecognitionException {
        // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:1:8: ( T__10 | T__11 | T__12 | T__13 | T__14 | T__15 | ROOTAGE_QUALIFIER | DECIMAL_NUMBER | UNQUOTED_ALPHA_TEXT | QUOTED_TEXT | NESTED_ML_COMMENT | WS )
        int alt5=12;
        int LA5_0 = input.LA(1);

        if ( (LA5_0=='#') ) {
            alt5=1;
        }
        else if ( (LA5_0=='(') ) {
            alt5=2;
        }
        else if ( (LA5_0==')') ) {
            alt5=3;
        }
        else if ( (LA5_0==',') ) {
            alt5=4;
        }
        else if ( (LA5_0==':') ) {
            alt5=5;
        }
        else if ( (LA5_0==';') ) {
            alt5=6;
        }
        else if ( (LA5_0=='[') ) {
            int LA5_7 = input.LA(2);

            if ( (LA5_7=='&') ) {
                int LA5_12 = input.LA(3);

                if ( (LA5_12=='R'||LA5_12=='U'||LA5_12=='r'||LA5_12=='u') ) {
                    int LA5_14 = input.LA(4);

                    if ( (LA5_14==']') ) {
                        alt5=7;
                    }
                    else if ( ((LA5_14 >= '\u0000' && LA5_14 <= '\\')||(LA5_14 >= '^' && LA5_14 <= '\uFFFF')) ) {
                        alt5=11;
                    }
                    else {
                        NoViableAltException nvae =
                            new NoViableAltException("", 5, 14, input);

                        throw nvae;

                    }
                }
                else if ( ((LA5_12 >= '\u0000' && LA5_12 <= 'Q')||(LA5_12 >= 'S' && LA5_12 <= 'T')||(LA5_12 >= 'V' && LA5_12 <= 'q')||(LA5_12 >= 's' && LA5_12 <= 't')||(LA5_12 >= 'v' && LA5_12 <= '\uFFFF')) ) {
                    alt5=11;
                }
                else {
                    NoViableAltException nvae =
                        new NoViableAltException("", 5, 12, input);

                    throw nvae;

                }
            }
            else if ( ((LA5_7 >= '\u0000' && LA5_7 <= '%')||(LA5_7 >= '\'' && LA5_7 <= '\uFFFF')) ) {
                alt5=11;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 5, 7, input);

                throw nvae;

            }
        }
        else if ( (LA5_0=='.'||(LA5_0 >= '0' && LA5_0 <= '9')) ) {
            alt5=8;
        }
        else if ( ((LA5_0 >= '\u0000' && LA5_0 <= '\b')||(LA5_0 >= '\u000B' && LA5_0 <= '\f')||(LA5_0 >= '\u000E' && LA5_0 <= '\u001F')||(LA5_0 >= '!' && LA5_0 <= '\"')||(LA5_0 >= '$' && LA5_0 <= '&')||(LA5_0 >= '*' && LA5_0 <= '+')||LA5_0=='-'||LA5_0=='/'||(LA5_0 >= '<' && LA5_0 <= 'Z')||LA5_0=='\\'||(LA5_0 >= '^' && LA5_0 <= '\uFFFF')) ) {
            alt5=9;
        }
        else if ( (LA5_0=='\'') ) {
            alt5=10;
        }
        else if ( ((LA5_0 >= '\t' && LA5_0 <= '\n')||LA5_0=='\r'||LA5_0==' ') ) {
            alt5=12;
        }
        else {
            NoViableAltException nvae =
                new NoViableAltException("", 5, 0, input);

            throw nvae;

        }
        switch (alt5) {
            case 1 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:1:10: T__10
                {
                mT__10(); 


                }
                break;
            case 2 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:1:16: T__11
                {
                mT__11(); 


                }
                break;
            case 3 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:1:22: T__12
                {
                mT__12(); 


                }
                break;
            case 4 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:1:28: T__13
                {
                mT__13(); 


                }
                break;
            case 5 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:1:34: T__14
                {
                mT__14(); 


                }
                break;
            case 6 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:1:40: T__15
                {
                mT__15(); 


                }
                break;
            case 7 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:1:46: ROOTAGE_QUALIFIER
                {
                mROOTAGE_QUALIFIER(); 


                }
                break;
            case 8 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:1:64: DECIMAL_NUMBER
                {
                mDECIMAL_NUMBER(); 


                }
                break;
            case 9 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:1:79: UNQUOTED_ALPHA_TEXT
                {
                mUNQUOTED_ALPHA_TEXT(); 


                }
                break;
            case 10 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:1:99: QUOTED_TEXT
                {
                mQUOTED_TEXT(); 


                }
                break;
            case 11 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:1:111: NESTED_ML_COMMENT
                {
                mNESTED_ML_COMMENT(); 


                }
                break;
            case 12 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:1:129: WS
                {
                mWS(); 


                }
                break;

        }

    }


 

}
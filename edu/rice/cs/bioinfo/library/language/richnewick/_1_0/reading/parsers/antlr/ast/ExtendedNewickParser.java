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

package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast;
import java.util.LinkedList;



import org.antlr.runtime.*;

import java.util.List;

@SuppressWarnings({"all", "warnings", "unchecked"})
public class ExtendedNewickParser extends Parser {
    public static final String[] tokenNames = new String[] {
        "<invalid>", "<EOR>", "<DOWN>", "<UP>", "DECIMAL_NUMBER", "NESTED_ML_COMMENT", "QUOTED_TEXT", "ROOTAGE_QUALIFIER", "UNQUOTED_ALPHA_TEXT", "WS", "'#'", "'('", "')'", "','", "':'", "';'"
    };

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
    public Parser[] getDelegates() {
        return new Parser[] {};
    }

    // delegators


    public ExtendedNewickParser(TokenStream input) {
        this(input, new RecognizerSharedState());
    }
    public ExtendedNewickParser(TokenStream input, RecognizerSharedState state) {
        super(input, state);
    }

    public String[] getTokenNames() { return ExtendedNewickParser.tokenNames; }
    public String getGrammarFileName() { return "D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g"; }



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

       
    public void displayRecognitionError(String[] tokenNames, RecognitionException e) 
    {
            errors.add(new ErrorWrapper(getErrorMessage(e, tokenNames), e.line, e.c));
    }
        
    public List<ErrorWrapper> getErrors() 
    {
            return errors;
    }

    ParseStack stack = new ParseStackAction();

    ParseStack getParseStack()
    {
    	return stack;
    }





    // $ANTLR start "networks"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:53:1: networks : ( network )* ;
    public final void networks() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:53:9: ( ( network )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:53:11: ( network )*
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:53:11: ( network )*
            loop1:
            do {
                int alt1=2;
                int LA1_0 = input.LA(1);

                if ( (LA1_0==DECIMAL_NUMBER||(LA1_0 >= QUOTED_TEXT && LA1_0 <= UNQUOTED_ALPHA_TEXT)||(LA1_0 >= 10 && LA1_0 <= 15)) ) {
                    alt1=1;
                }


                switch (alt1) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:53:11: network
            	    {
            	    pushFollow(FOLLOW_network_in_networks33);
            	    network();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop1;
                }
            } while (true);


             stack.pushNetworks(); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return ;
    }
    // $ANTLR end "networks"



    // $ANTLR start "network"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:55:1: network : (rq= ROOTAGE_QUALIFIER )? (dl= descendant_list )? network_info ';' ;
    public final void network() throws RecognitionException {
        Token rq=null;
        ExtendedNewickParser.descendant_list_return dl =null;


        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:55:9: ( (rq= ROOTAGE_QUALIFIER )? (dl= descendant_list )? network_info ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:55:11: (rq= ROOTAGE_QUALIFIER )? (dl= descendant_list )? network_info ';'
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:55:13: (rq= ROOTAGE_QUALIFIER )?
            int alt2=2;
            int LA2_0 = input.LA(1);

            if ( (LA2_0==ROOTAGE_QUALIFIER) ) {
                alt2=1;
            }
            switch (alt2) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:55:13: rq= ROOTAGE_QUALIFIER
                    {
                    rq=(Token)match(input,ROOTAGE_QUALIFIER,FOLLOW_ROOTAGE_QUALIFIER_in_network46); 

                    }
                    break;

            }


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:55:35: (dl= descendant_list )?
            int alt3=2;
            int LA3_0 = input.LA(1);

            if ( (LA3_0==11) ) {
                alt3=1;
            }
            switch (alt3) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:55:35: dl= descendant_list
                    {
                    pushFollow(FOLLOW_descendant_list_in_network51);
                    dl=descendant_list();

                    state._fsp--;


                    }
                    break;

            }


            pushFollow(FOLLOW_network_info_in_network54);
            network_info();

            state._fsp--;


            match(input,15,FOLLOW_15_in_network56); 

             stack.pushNetwork(rq, (dl!=null?((Token)dl.start):null) != null); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return ;
    }
    // $ANTLR end "network"


    public static class descendant_list_return extends ParserRuleReturnScope {
    };


    // $ANTLR start "descendant_list"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:57:1: descendant_list : '(' subtree ( ',' subtree )* ')' ;
    public final ExtendedNewickParser.descendant_list_return descendant_list() throws RecognitionException {
        ExtendedNewickParser.descendant_list_return retval = new ExtendedNewickParser.descendant_list_return();
        retval.start = input.LT(1);


         int subTreeCount = 1; 
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:57:48: ( '(' subtree ( ',' subtree )* ')' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:57:50: '(' subtree ( ',' subtree )* ')'
            {
            match(input,11,FOLLOW_11_in_descendant_list72); 

            pushFollow(FOLLOW_subtree_in_descendant_list74);
            subtree();

            state._fsp--;


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:57:62: ( ',' subtree )*
            loop4:
            do {
                int alt4=2;
                int LA4_0 = input.LA(1);

                if ( (LA4_0==13) ) {
                    alt4=1;
                }


                switch (alt4) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:57:63: ',' subtree
            	    {
            	    match(input,13,FOLLOW_13_in_descendant_list77); 

            	    pushFollow(FOLLOW_subtree_in_descendant_list79);
            	    subtree();

            	    state._fsp--;


            	    subTreeCount++;

            	    }
            	    break;

            	default :
            	    break loop4;
                }
            } while (true);


            match(input,12,FOLLOW_12_in_descendant_list86); 

             stack.pushDescendantList(subTreeCount); 

            }

            retval.stop = input.LT(-1);


        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return retval;
    }
    // $ANTLR end "descendant_list"



    // $ANTLR start "subtree"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:59:1: subtree : ( descendant_list network_info | network_info );
    public final void subtree() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:59:9: ( descendant_list network_info | network_info )
            int alt5=2;
            int LA5_0 = input.LA(1);

            if ( (LA5_0==11) ) {
                alt5=1;
            }
            else if ( (LA5_0==DECIMAL_NUMBER||LA5_0==QUOTED_TEXT||LA5_0==UNQUOTED_ALPHA_TEXT||LA5_0==10||(LA5_0 >= 12 && LA5_0 <= 15)) ) {
                alt5=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 5, 0, input);

                throw nvae;

            }
            switch (alt5) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:59:11: descendant_list network_info
                    {
                    pushFollow(FOLLOW_descendant_list_in_subtree96);
                    descendant_list();

                    state._fsp--;


                    pushFollow(FOLLOW_network_info_in_subtree98);
                    network_info();

                    state._fsp--;


                     stack.pushSubtree(true); 

                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:60:11: network_info
                    {
                    pushFollow(FOLLOW_network_info_in_subtree112);
                    network_info();

                    state._fsp--;


                     stack.pushSubtree(false); 

                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return ;
    }
    // $ANTLR end "subtree"



    // $ANTLR start "network_info"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:62:1: network_info : ( (nl= node_label )? (hn= hybrid_node_qualifier )? (bl= branch_length )? | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length support | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length support probability | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length ':' probability | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' support (pb= probability )? | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' ':' probability );
    public final void network_info() throws RecognitionException {
        ExtendedNewickParser.node_label_return nl =null;

        ExtendedNewickParser.hybrid_node_qualifier_return hn =null;

        ExtendedNewickParser.branch_length_return bl =null;

        ExtendedNewickParser.probability_return pb =null;


        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:63:2: ( (nl= node_label )? (hn= hybrid_node_qualifier )? (bl= branch_length )? | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length support | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length support probability | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length ':' probability | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' support (pb= probability )? | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' ':' probability )
            int alt20=6;
            alt20 = dfa20.predict(input);
            switch (alt20) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:63:4: (nl= node_label )? (hn= hybrid_node_qualifier )? (bl= branch_length )?
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:63:6: (nl= node_label )?
                    int alt6=2;
                    int LA6_0 = input.LA(1);

                    if ( (LA6_0==DECIMAL_NUMBER||LA6_0==QUOTED_TEXT||LA6_0==UNQUOTED_ALPHA_TEXT) ) {
                        alt6=1;
                    }
                    switch (alt6) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:63:6: nl= node_label
                            {
                            pushFollow(FOLLOW_node_label_in_network_info141);
                            nl=node_label();

                            state._fsp--;


                            }
                            break;

                    }


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:63:21: (hn= hybrid_node_qualifier )?
                    int alt7=2;
                    int LA7_0 = input.LA(1);

                    if ( (LA7_0==10) ) {
                        alt7=1;
                    }
                    switch (alt7) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:63:21: hn= hybrid_node_qualifier
                            {
                            pushFollow(FOLLOW_hybrid_node_qualifier_in_network_info146);
                            hn=hybrid_node_qualifier();

                            state._fsp--;


                            }
                            break;

                    }


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:63:47: (bl= branch_length )?
                    int alt8=2;
                    int LA8_0 = input.LA(1);

                    if ( (LA8_0==14) ) {
                        alt8=1;
                    }
                    switch (alt8) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:63:47: bl= branch_length
                            {
                            pushFollow(FOLLOW_branch_length_in_network_info151);
                            bl=branch_length();

                            state._fsp--;


                            }
                            break;

                    }


                     stack.pushNetworkInfo((nl!=null?((Token)nl.start):null) != null, (hn!=null?((Token)hn.start):null) != null, (bl!=null?((Token)bl.start):null) != null, false, false); 

                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:64:4: (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length support
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:64:6: (nl= node_label )?
                    int alt9=2;
                    int LA9_0 = input.LA(1);

                    if ( (LA9_0==DECIMAL_NUMBER||LA9_0==QUOTED_TEXT||LA9_0==UNQUOTED_ALPHA_TEXT) ) {
                        alt9=1;
                    }
                    switch (alt9) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:64:6: nl= node_label
                            {
                            pushFollow(FOLLOW_node_label_in_network_info180);
                            nl=node_label();

                            state._fsp--;


                            }
                            break;

                    }


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:64:21: (hn= hybrid_node_qualifier )?
                    int alt10=2;
                    int LA10_0 = input.LA(1);

                    if ( (LA10_0==10) ) {
                        alt10=1;
                    }
                    switch (alt10) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:64:21: hn= hybrid_node_qualifier
                            {
                            pushFollow(FOLLOW_hybrid_node_qualifier_in_network_info185);
                            hn=hybrid_node_qualifier();

                            state._fsp--;


                            }
                            break;

                    }


                    pushFollow(FOLLOW_branch_length_in_network_info191);
                    branch_length();

                    state._fsp--;


                    pushFollow(FOLLOW_support_in_network_info194);
                    support();

                    state._fsp--;


                     stack.pushNetworkInfo((nl!=null?((Token)nl.start):null) != null, (hn!=null?((Token)hn.start):null) != null, true, true, false); 

                    }
                    break;
                case 3 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:65:4: (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length support probability
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:65:6: (nl= node_label )?
                    int alt11=2;
                    int LA11_0 = input.LA(1);

                    if ( (LA11_0==DECIMAL_NUMBER||LA11_0==QUOTED_TEXT||LA11_0==UNQUOTED_ALPHA_TEXT) ) {
                        alt11=1;
                    }
                    switch (alt11) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:65:6: nl= node_label
                            {
                            pushFollow(FOLLOW_node_label_in_network_info212);
                            nl=node_label();

                            state._fsp--;


                            }
                            break;

                    }


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:65:21: (hn= hybrid_node_qualifier )?
                    int alt12=2;
                    int LA12_0 = input.LA(1);

                    if ( (LA12_0==10) ) {
                        alt12=1;
                    }
                    switch (alt12) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:65:21: hn= hybrid_node_qualifier
                            {
                            pushFollow(FOLLOW_hybrid_node_qualifier_in_network_info217);
                            hn=hybrid_node_qualifier();

                            state._fsp--;


                            }
                            break;

                    }


                    pushFollow(FOLLOW_branch_length_in_network_info220);
                    branch_length();

                    state._fsp--;


                    pushFollow(FOLLOW_support_in_network_info223);
                    support();

                    state._fsp--;


                    pushFollow(FOLLOW_probability_in_network_info225);
                    probability();

                    state._fsp--;


                     stack.pushNetworkInfo((nl!=null?((Token)nl.start):null) != null, (hn!=null?((Token)hn.start):null) != null, true, true, true); 

                    }
                    break;
                case 4 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:66:4: (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length ':' probability
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:66:6: (nl= node_label )?
                    int alt13=2;
                    int LA13_0 = input.LA(1);

                    if ( (LA13_0==DECIMAL_NUMBER||LA13_0==QUOTED_TEXT||LA13_0==UNQUOTED_ALPHA_TEXT) ) {
                        alt13=1;
                    }
                    switch (alt13) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:66:6: nl= node_label
                            {
                            pushFollow(FOLLOW_node_label_in_network_info234);
                            nl=node_label();

                            state._fsp--;


                            }
                            break;

                    }


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:66:21: (hn= hybrid_node_qualifier )?
                    int alt14=2;
                    int LA14_0 = input.LA(1);

                    if ( (LA14_0==10) ) {
                        alt14=1;
                    }
                    switch (alt14) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:66:21: hn= hybrid_node_qualifier
                            {
                            pushFollow(FOLLOW_hybrid_node_qualifier_in_network_info239);
                            hn=hybrid_node_qualifier();

                            state._fsp--;


                            }
                            break;

                    }


                    pushFollow(FOLLOW_branch_length_in_network_info242);
                    branch_length();

                    state._fsp--;


                    match(input,14,FOLLOW_14_in_network_info245); 

                    pushFollow(FOLLOW_probability_in_network_info247);
                    probability();

                    state._fsp--;


                     stack.pushNetworkInfo((nl!=null?((Token)nl.start):null) != null, (hn!=null?((Token)hn.start):null) != null, true, false, true); 

                    }
                    break;
                case 5 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:67:4: (nl= node_label )? (hn= hybrid_node_qualifier )? ':' support (pb= probability )?
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:67:6: (nl= node_label )?
                    int alt15=2;
                    int LA15_0 = input.LA(1);

                    if ( (LA15_0==DECIMAL_NUMBER||LA15_0==QUOTED_TEXT||LA15_0==UNQUOTED_ALPHA_TEXT) ) {
                        alt15=1;
                    }
                    switch (alt15) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:67:6: nl= node_label
                            {
                            pushFollow(FOLLOW_node_label_in_network_info262);
                            nl=node_label();

                            state._fsp--;


                            }
                            break;

                    }


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:67:21: (hn= hybrid_node_qualifier )?
                    int alt16=2;
                    int LA16_0 = input.LA(1);

                    if ( (LA16_0==10) ) {
                        alt16=1;
                    }
                    switch (alt16) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:67:21: hn= hybrid_node_qualifier
                            {
                            pushFollow(FOLLOW_hybrid_node_qualifier_in_network_info267);
                            hn=hybrid_node_qualifier();

                            state._fsp--;


                            }
                            break;

                    }


                    match(input,14,FOLLOW_14_in_network_info270); 

                    pushFollow(FOLLOW_support_in_network_info273);
                    support();

                    state._fsp--;


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:67:60: (pb= probability )?
                    int alt17=2;
                    int LA17_0 = input.LA(1);

                    if ( (LA17_0==14) ) {
                        alt17=1;
                    }
                    switch (alt17) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:67:60: pb= probability
                            {
                            pushFollow(FOLLOW_probability_in_network_info277);
                            pb=probability();

                            state._fsp--;


                            }
                            break;

                    }


                     stack.pushNetworkInfo((nl!=null?((Token)nl.start):null) != null, (hn!=null?((Token)hn.start):null) != null, false, true, (pb!=null?((Token)pb.start):null) != null); 

                    }
                    break;
                case 6 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:68:5: (nl= node_label )? (hn= hybrid_node_qualifier )? ':' ':' probability
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:68:7: (nl= node_label )?
                    int alt18=2;
                    int LA18_0 = input.LA(1);

                    if ( (LA18_0==DECIMAL_NUMBER||LA18_0==QUOTED_TEXT||LA18_0==UNQUOTED_ALPHA_TEXT) ) {
                        alt18=1;
                    }
                    switch (alt18) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:68:7: nl= node_label
                            {
                            pushFollow(FOLLOW_node_label_in_network_info294);
                            nl=node_label();

                            state._fsp--;


                            }
                            break;

                    }


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:68:22: (hn= hybrid_node_qualifier )?
                    int alt19=2;
                    int LA19_0 = input.LA(1);

                    if ( (LA19_0==10) ) {
                        alt19=1;
                    }
                    switch (alt19) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:68:22: hn= hybrid_node_qualifier
                            {
                            pushFollow(FOLLOW_hybrid_node_qualifier_in_network_info299);
                            hn=hybrid_node_qualifier();

                            state._fsp--;


                            }
                            break;

                    }


                    match(input,14,FOLLOW_14_in_network_info302); 

                    match(input,14,FOLLOW_14_in_network_info305); 

                    pushFollow(FOLLOW_probability_in_network_info307);
                    probability();

                    state._fsp--;


                     stack.pushNetworkInfo((nl!=null?((Token)nl.start):null) != null, (hn!=null?((Token)hn.start):null) != null, false, false, true); 

                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return ;
    }
    // $ANTLR end "network_info"


    public static class branch_length_return extends ParserRuleReturnScope {
    };


    // $ANTLR start "branch_length"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:70:1: branch_length : edge_label ;
    public final ExtendedNewickParser.branch_length_return branch_length() throws RecognitionException {
        ExtendedNewickParser.branch_length_return retval = new ExtendedNewickParser.branch_length_return();
        retval.start = input.LT(1);


        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:71:2: ( edge_label )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:71:4: edge_label
            {
            pushFollow(FOLLOW_edge_label_in_branch_length328);
            edge_label();

            state._fsp--;


             stack.pushBranchLength(); 

            }

            retval.stop = input.LT(-1);


        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return retval;
    }
    // $ANTLR end "branch_length"



    // $ANTLR start "support"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:73:1: support : edge_label ;
    public final void support() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:74:2: ( edge_label )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:74:4: edge_label
            {
            pushFollow(FOLLOW_edge_label_in_support340);
            edge_label();

            state._fsp--;


             stack.pushSupport(); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return ;
    }
    // $ANTLR end "support"


    public static class probability_return extends ParserRuleReturnScope {
    };


    // $ANTLR start "probability"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:76:1: probability : edge_label ;
    public final ExtendedNewickParser.probability_return probability() throws RecognitionException {
        ExtendedNewickParser.probability_return retval = new ExtendedNewickParser.probability_return();
        retval.start = input.LT(1);


        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:77:2: ( edge_label )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:77:4: edge_label
            {
            pushFollow(FOLLOW_edge_label_in_probability352);
            edge_label();

            state._fsp--;


             stack.pushProbability(); 

            }

            retval.stop = input.LT(-1);


        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return retval;
    }
    // $ANTLR end "probability"


    public static class node_label_return extends ParserRuleReturnScope {
    };


    // $ANTLR start "node_label"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:79:1: node_label : text ;
    public final ExtendedNewickParser.node_label_return node_label() throws RecognitionException {
        ExtendedNewickParser.node_label_return retval = new ExtendedNewickParser.node_label_return();
        retval.start = input.LT(1);


        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:80:2: ( text )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:80:4: text
            {
            pushFollow(FOLLOW_text_in_node_label363);
            text();

            state._fsp--;


             stack.pushNodeLabel(); 

            }

            retval.stop = input.LT(-1);


        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return retval;
    }
    // $ANTLR end "node_label"


    public static class hybrid_node_qualifier_return extends ParserRuleReturnScope {
    };


    // $ANTLR start "hybrid_node_qualifier"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:82:1: hybrid_node_qualifier : '#' (type= UNQUOTED_ALPHA_TEXT )? hybridNodeIndex= DECIMAL_NUMBER ;
    public final ExtendedNewickParser.hybrid_node_qualifier_return hybrid_node_qualifier() throws RecognitionException {
        ExtendedNewickParser.hybrid_node_qualifier_return retval = new ExtendedNewickParser.hybrid_node_qualifier_return();
        retval.start = input.LT(1);


        Token type=null;
        Token hybridNodeIndex=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:83:2: ( '#' (type= UNQUOTED_ALPHA_TEXT )? hybridNodeIndex= DECIMAL_NUMBER )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:83:4: '#' (type= UNQUOTED_ALPHA_TEXT )? hybridNodeIndex= DECIMAL_NUMBER
            {
            match(input,10,FOLLOW_10_in_hybrid_node_qualifier376); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:83:12: (type= UNQUOTED_ALPHA_TEXT )?
            int alt21=2;
            int LA21_0 = input.LA(1);

            if ( (LA21_0==UNQUOTED_ALPHA_TEXT) ) {
                alt21=1;
            }
            switch (alt21) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:83:12: type= UNQUOTED_ALPHA_TEXT
                    {
                    type=(Token)match(input,UNQUOTED_ALPHA_TEXT,FOLLOW_UNQUOTED_ALPHA_TEXT_in_hybrid_node_qualifier380); 

                    }
                    break;

            }


            hybridNodeIndex=(Token)match(input,DECIMAL_NUMBER,FOLLOW_DECIMAL_NUMBER_in_hybrid_node_qualifier385); 

             stack.pushUnquotedText((hybridNodeIndex!=null?hybridNodeIndex.getText():null), (hybridNodeIndex!=null?hybridNodeIndex.getLine():0), (hybridNodeIndex!=null?hybridNodeIndex.getCharPositionInLine():0)); 
            		  if(type!=null)stack.pushUnquotedText((type!=null?type.getText():null), (type!=null?type.getLine():0), (type!=null?type.getCharPositionInLine():0)); stack.pushHybridNodeQualifier(type, hybridNodeIndex); 

            }

            retval.stop = input.LT(-1);


        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return retval;
    }
    // $ANTLR end "hybrid_node_qualifier"



    // $ANTLR start "edge_label"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:87:1: edge_label : ':' text ;
    public final void edge_label() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:88:2: ( ':' text )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:88:4: ':' text
            {
            match(input,14,FOLLOW_14_in_edge_label399); 

            pushFollow(FOLLOW_text_in_edge_label401);
            text();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return ;
    }
    // $ANTLR end "edge_label"



    // $ANTLR start "text"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:91:1: text : (q= QUOTED_TEXT | (t= ( UNQUOTED_ALPHA_TEXT | DECIMAL_NUMBER ) )+ );
    public final void text() throws RecognitionException {
        Token q=null;
        Token t=null;

         String str = ""; int ln = -1; int cn = -1; 
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:91:58: (q= QUOTED_TEXT | (t= ( UNQUOTED_ALPHA_TEXT | DECIMAL_NUMBER ) )+ )
            int alt23=2;
            int LA23_0 = input.LA(1);

            if ( (LA23_0==QUOTED_TEXT) ) {
                alt23=1;
            }
            else if ( (LA23_0==DECIMAL_NUMBER||LA23_0==UNQUOTED_ALPHA_TEXT) ) {
                alt23=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 23, 0, input);

                throw nvae;

            }
            switch (alt23) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:92:14: q= QUOTED_TEXT
                    {
                    q=(Token)match(input,QUOTED_TEXT,FOLLOW_QUOTED_TEXT_in_text431); 

                     stack.pushQuotedText(q); 

                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:93:7: (t= ( UNQUOTED_ALPHA_TEXT | DECIMAL_NUMBER ) )+
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:93:7: (t= ( UNQUOTED_ALPHA_TEXT | DECIMAL_NUMBER ) )+
                    int cnt22=0;
                    loop22:
                    do {
                        int alt22=2;
                        int LA22_0 = input.LA(1);

                        if ( (LA22_0==DECIMAL_NUMBER||LA22_0==UNQUOTED_ALPHA_TEXT) ) {
                            alt22=1;
                        }


                        switch (alt22) {
                    	case 1 :
                    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\ExtendedNewick.g:93:8: t= ( UNQUOTED_ALPHA_TEXT | DECIMAL_NUMBER )
                    	    {
                    	    t=(Token)input.LT(1);

                    	    if ( input.LA(1)==DECIMAL_NUMBER||input.LA(1)==UNQUOTED_ALPHA_TEXT ) {
                    	        input.consume();
                    	        state.errorRecovery=false;
                    	    }
                    	    else {
                    	        MismatchedSetException mse = new MismatchedSetException(null,input);
                    	        throw mse;
                    	    }


                    	    str+=(t!=null?t.getText():null); ln = ln == -1 ? t.getLine() : ln; cn = cn == -1 ? t.getCharPositionInLine() : cn;

                    	    }
                    	    break;

                    	default :
                    	    if ( cnt22 >= 1 ) break loop22;
                                EarlyExitException eee =
                                    new EarlyExitException(22, input);
                                throw eee;
                        }
                        cnt22++;
                    } while (true);


                     
                    				  

                            				str = str.trim();

                            				if(str.contains(" ") || str.contains("\n") || str.contains("\r") || str.contains("\t") )
                            				{
                                					errors.add(new ErrorWrapper("Invalid whitespace in node or edge label '" + str + "'.", ln, cn));
                            				}
                    				  
                    				  	stack.pushUnquotedText(str, ln, cn); 
                    				  

                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return ;
    }
    // $ANTLR end "text"

    // Delegated rules


    protected DFA20 dfa20 = new DFA20(this);
    static final String DFA20_eotS =
        "\23\uffff";
    static final String DFA20_eofS =
        "\23\uffff";
    static final String DFA20_minS =
        "\1\4\1\12\3\4\1\uffff\1\4\1\14\1\4\1\14\1\4\2\uffff\1\4\1\14\1\4"+
        "\3\uffff";
    static final String DFA20_maxS =
        "\3\17\1\10\1\16\1\uffff\1\4\1\17\1\16\2\17\2\uffff\1\16\2\17\3\uffff";
    static final String DFA20_acceptS =
        "\5\uffff\1\1\5\uffff\1\6\1\5\3\uffff\1\4\1\2\1\3";
    static final String DFA20_specialS =
        "\23\uffff}>";
    static final String[] DFA20_transitionS = {
            "\1\2\1\uffff\1\1\1\uffff\1\2\1\uffff\1\3\1\uffff\2\5\1\4\1\5",
            "\1\3\1\uffff\2\5\1\4\1\5",
            "\1\2\3\uffff\1\2\1\uffff\1\3\1\uffff\2\5\1\4\1\5",
            "\1\7\3\uffff\1\6",
            "\1\12\1\uffff\1\11\1\uffff\1\12\5\uffff\1\10",
            "",
            "\1\7",
            "\2\5\1\4\1\5",
            "\1\14\1\uffff\1\14\1\uffff\1\14\5\uffff\1\13",
            "\2\5\1\15\1\5",
            "\1\12\3\uffff\1\12\3\uffff\2\5\1\15\1\5",
            "",
            "",
            "\1\17\1\uffff\1\16\1\uffff\1\17\5\uffff\1\20",
            "\2\21\1\22\1\21",
            "\1\17\3\uffff\1\17\3\uffff\2\21\1\22\1\21",
            "",
            "",
            ""
    };

    static final short[] DFA20_eot = DFA.unpackEncodedString(DFA20_eotS);
    static final short[] DFA20_eof = DFA.unpackEncodedString(DFA20_eofS);
    static final char[] DFA20_min = DFA.unpackEncodedStringToUnsignedChars(DFA20_minS);
    static final char[] DFA20_max = DFA.unpackEncodedStringToUnsignedChars(DFA20_maxS);
    static final short[] DFA20_accept = DFA.unpackEncodedString(DFA20_acceptS);
    static final short[] DFA20_special = DFA.unpackEncodedString(DFA20_specialS);
    static final short[][] DFA20_transition;

    static {
        int numStates = DFA20_transitionS.length;
        DFA20_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA20_transition[i] = DFA.unpackEncodedString(DFA20_transitionS[i]);
        }
    }

    class DFA20 extends DFA {

        public DFA20(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 20;
            this.eot = DFA20_eot;
            this.eof = DFA20_eof;
            this.min = DFA20_min;
            this.max = DFA20_max;
            this.accept = DFA20_accept;
            this.special = DFA20_special;
            this.transition = DFA20_transition;
        }
        public String getDescription() {
            return "62:1: network_info : ( (nl= node_label )? (hn= hybrid_node_qualifier )? (bl= branch_length )? | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length support | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length support probability | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length ':' probability | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' support (pb= probability )? | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' ':' probability );";
        }
    }
 

    public static final BitSet FOLLOW_network_in_networks33 = new BitSet(new long[]{0x000000000000CDD2L});
    public static final BitSet FOLLOW_ROOTAGE_QUALIFIER_in_network46 = new BitSet(new long[]{0x000000000000CD50L});
    public static final BitSet FOLLOW_descendant_list_in_network51 = new BitSet(new long[]{0x000000000000C550L});
    public static final BitSet FOLLOW_network_info_in_network54 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_15_in_network56 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_11_in_descendant_list72 = new BitSet(new long[]{0x0000000000004D50L});
    public static final BitSet FOLLOW_subtree_in_descendant_list74 = new BitSet(new long[]{0x0000000000003000L});
    public static final BitSet FOLLOW_13_in_descendant_list77 = new BitSet(new long[]{0x0000000000004D50L});
    public static final BitSet FOLLOW_subtree_in_descendant_list79 = new BitSet(new long[]{0x0000000000003000L});
    public static final BitSet FOLLOW_12_in_descendant_list86 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_descendant_list_in_subtree96 = new BitSet(new long[]{0x0000000000004550L});
    public static final BitSet FOLLOW_network_info_in_subtree98 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_network_info_in_subtree112 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_node_label_in_network_info141 = new BitSet(new long[]{0x0000000000004402L});
    public static final BitSet FOLLOW_hybrid_node_qualifier_in_network_info146 = new BitSet(new long[]{0x0000000000004002L});
    public static final BitSet FOLLOW_branch_length_in_network_info151 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_node_label_in_network_info180 = new BitSet(new long[]{0x0000000000004400L});
    public static final BitSet FOLLOW_hybrid_node_qualifier_in_network_info185 = new BitSet(new long[]{0x0000000000004000L});
    public static final BitSet FOLLOW_branch_length_in_network_info191 = new BitSet(new long[]{0x0000000000004000L});
    public static final BitSet FOLLOW_support_in_network_info194 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_node_label_in_network_info212 = new BitSet(new long[]{0x0000000000004400L});
    public static final BitSet FOLLOW_hybrid_node_qualifier_in_network_info217 = new BitSet(new long[]{0x0000000000004000L});
    public static final BitSet FOLLOW_branch_length_in_network_info220 = new BitSet(new long[]{0x0000000000004000L});
    public static final BitSet FOLLOW_support_in_network_info223 = new BitSet(new long[]{0x0000000000004000L});
    public static final BitSet FOLLOW_probability_in_network_info225 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_node_label_in_network_info234 = new BitSet(new long[]{0x0000000000004400L});
    public static final BitSet FOLLOW_hybrid_node_qualifier_in_network_info239 = new BitSet(new long[]{0x0000000000004000L});
    public static final BitSet FOLLOW_branch_length_in_network_info242 = new BitSet(new long[]{0x0000000000004000L});
    public static final BitSet FOLLOW_14_in_network_info245 = new BitSet(new long[]{0x0000000000004000L});
    public static final BitSet FOLLOW_probability_in_network_info247 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_node_label_in_network_info262 = new BitSet(new long[]{0x0000000000004400L});
    public static final BitSet FOLLOW_hybrid_node_qualifier_in_network_info267 = new BitSet(new long[]{0x0000000000004000L});
    public static final BitSet FOLLOW_14_in_network_info270 = new BitSet(new long[]{0x0000000000004000L});
    public static final BitSet FOLLOW_support_in_network_info273 = new BitSet(new long[]{0x0000000000004002L});
    public static final BitSet FOLLOW_probability_in_network_info277 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_node_label_in_network_info294 = new BitSet(new long[]{0x0000000000004400L});
    public static final BitSet FOLLOW_hybrid_node_qualifier_in_network_info299 = new BitSet(new long[]{0x0000000000004000L});
    public static final BitSet FOLLOW_14_in_network_info302 = new BitSet(new long[]{0x0000000000004000L});
    public static final BitSet FOLLOW_14_in_network_info305 = new BitSet(new long[]{0x0000000000004000L});
    public static final BitSet FOLLOW_probability_in_network_info307 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_edge_label_in_branch_length328 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_edge_label_in_support340 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_edge_label_in_probability352 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_text_in_node_label363 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_10_in_hybrid_node_qualifier376 = new BitSet(new long[]{0x0000000000000110L});
    public static final BitSet FOLLOW_UNQUOTED_ALPHA_TEXT_in_hybrid_node_qualifier380 = new BitSet(new long[]{0x0000000000000010L});
    public static final BitSet FOLLOW_DECIMAL_NUMBER_in_hybrid_node_qualifier385 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_14_in_edge_label399 = new BitSet(new long[]{0x0000000000000150L});
    public static final BitSet FOLLOW_text_in_edge_label401 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_QUOTED_TEXT_in_text431 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_set_in_text447 = new BitSet(new long[]{0x0000000000000112L});

}
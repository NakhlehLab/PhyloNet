// $ANTLR 3.4 D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g 2011-08-31 18:58:15

package edu.rice.bioinfo.library.language.pyson._1_0.parsers.antlr.ast;
import java.util.LinkedList;


import org.antlr.runtime.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;

@SuppressWarnings({"all", "warnings", "unchecked"})
public class PySONParser extends Parser {
    public static final String[] tokenNames = new String[] {
        "<invalid>", "<EOR>", "<DOWN>", "<UP>", "BEGIN", "DEFAULT_INDICATOR", "END", "ID", "NESTED_ML_COMMENT", "NETWORK", "NETWORKS", "START", "TRANSLATE", "TREE", "TREES", "UTREE", "WS", "','", "';'", "'='"
    };

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
    public Parser[] getDelegates() {
        return new Parser[] {};
    }

    // delegators


    public PySONParser(TokenStream input) {
        this(input, new RecognizerSharedState());
    }
    public PySONParser(TokenStream input, RecognizerSharedState state) {
        super(input, state);
    }

    public String[] getTokenNames() { return PySONParser.tokenNames; }
    public String getGrammarFileName() { return "D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g"; }



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





    // $ANTLR start "blocks"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:54:1: blocks : START ( block )* ;
    public final void blocks() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:54:9: ( START ( block )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:54:11: START ( block )*
            {
            match(input,START,FOLLOW_START_in_blocks36); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:54:17: ( block )*
            loop1:
            do {
                int alt1=2;
                int LA1_0 = input.LA(1);

                if ( (LA1_0==BEGIN) ) {
                    alt1=1;
                }


                switch (alt1) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:54:17: block
            	    {
            	    pushFollow(FOLLOW_block_in_blocks38);
            	    block();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop1;
                }
            } while (true);


             stack.pushBlocks(); 

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
    // $ANTLR end "blocks"



    // $ANTLR start "block"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:56:1: block : ( BEGIN networks_block_body END | BEGIN trees_block_body END | BEGIN non_networks_block_body END );
    public final void block() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:56:7: ( BEGIN networks_block_body END | BEGIN trees_block_body END | BEGIN non_networks_block_body END )
            int alt2=3;
            int LA2_0 = input.LA(1);

            if ( (LA2_0==BEGIN) ) {
                switch ( input.LA(2) ) {
                case NETWORKS:
                    {
                    alt2=1;
                    }
                    break;
                case TREES:
                    {
                    alt2=2;
                    }
                    break;
                case BEGIN:
                case DEFAULT_INDICATOR:
                case END:
                case ID:
                case NESTED_ML_COMMENT:
                case NETWORK:
                case START:
                case TRANSLATE:
                case TREE:
                case UTREE:
                case WS:
                case 17:
                case 18:
                case 19:
                    {
                    alt2=3;
                    }
                    break;
                default:
                    NoViableAltException nvae =
                        new NoViableAltException("", 2, 1, input);

                    throw nvae;

                }

            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 2, 0, input);

                throw nvae;

            }
            switch (alt2) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:56:9: BEGIN networks_block_body END
                    {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block49); 

                    pushFollow(FOLLOW_networks_block_body_in_block51);
                    networks_block_body();

                    state._fsp--;


                    match(input,END,FOLLOW_END_in_block57); 

                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:57:4: BEGIN trees_block_body END
                    {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block62); 

                    pushFollow(FOLLOW_trees_block_body_in_block64);
                    trees_block_body();

                    state._fsp--;


                    match(input,END,FOLLOW_END_in_block73); 

                    }
                    break;
                case 3 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:58:4: BEGIN non_networks_block_body END
                    {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block78); 

                    pushFollow(FOLLOW_non_networks_block_body_in_block80);
                    non_networks_block_body();

                    state._fsp--;


                    match(input,END,FOLLOW_END_in_block82); 

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
    // $ANTLR end "block"



    // $ANTLR start "networks_block_body"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:60:1: networks_block_body : NETWORKS translation NETWORK rich_newick_assignment ;
    public final void networks_block_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:61:2: ( NETWORKS translation NETWORK rich_newick_assignment )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:62:3: NETWORKS translation NETWORK rich_newick_assignment
            {
            match(input,NETWORKS,FOLLOW_NETWORKS_in_networks_block_body93); 

            pushFollow(FOLLOW_translation_in_networks_block_body95);
            translation();

            state._fsp--;


            match(input,NETWORK,FOLLOW_NETWORK_in_networks_block_body97); 

            pushFollow(FOLLOW_rich_newick_assignment_in_networks_block_body99);
            rich_newick_assignment();

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
    // $ANTLR end "networks_block_body"



    // $ANTLR start "trees_block_body"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:64:1: trees_block_body : ( TREES translation ( tree_assigment )* | TREES ( tree_assigment )* );
    public final void trees_block_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:65:2: ( TREES translation ( tree_assigment )* | TREES ( tree_assigment )* )
            int alt5=2;
            int LA5_0 = input.LA(1);

            if ( (LA5_0==TREES) ) {
                int LA5_1 = input.LA(2);

                if ( (LA5_1==TRANSLATE) ) {
                    alt5=1;
                }
                else if ( (LA5_1==END||LA5_1==TREE||LA5_1==UTREE) ) {
                    alt5=2;
                }
                else {
                    NoViableAltException nvae =
                        new NoViableAltException("", 5, 1, input);

                    throw nvae;

                }
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 5, 0, input);

                throw nvae;

            }
            switch (alt5) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:65:4: TREES translation ( tree_assigment )*
                    {
                    match(input,TREES,FOLLOW_TREES_in_trees_block_body110); 

                    pushFollow(FOLLOW_translation_in_trees_block_body112);
                    translation();

                    state._fsp--;


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:65:22: ( tree_assigment )*
                    loop3:
                    do {
                        int alt3=2;
                        int LA3_0 = input.LA(1);

                        if ( (LA3_0==TREE||LA3_0==UTREE) ) {
                            alt3=1;
                        }


                        switch (alt3) {
                    	case 1 :
                    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:65:22: tree_assigment
                    	    {
                    	    pushFollow(FOLLOW_tree_assigment_in_trees_block_body114);
                    	    tree_assigment();

                    	    state._fsp--;


                    	    }
                    	    break;

                    	default :
                    	    break loop3;
                        }
                    } while (true);


                     stack.pushTreesBlockBody(true); 

                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:66:4: TREES ( tree_assigment )*
                    {
                    match(input,TREES,FOLLOW_TREES_in_trees_block_body122); 

                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:66:22: ( tree_assigment )*
                    loop4:
                    do {
                        int alt4=2;
                        int LA4_0 = input.LA(1);

                        if ( (LA4_0==TREE||LA4_0==UTREE) ) {
                            alt4=1;
                        }


                        switch (alt4) {
                    	case 1 :
                    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:66:22: tree_assigment
                    	    {
                    	    pushFollow(FOLLOW_tree_assigment_in_trees_block_body136);
                    	    tree_assigment();

                    	    state._fsp--;


                    	    }
                    	    break;

                    	default :
                    	    break loop4;
                        }
                    } while (true);


                     stack.pushTreesBlockBody(false); 

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
    // $ANTLR end "trees_block_body"



    // $ANTLR start "tree_assigment"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:68:1: tree_assigment : tr= ( TREE | UTREE ) rich_newick_assignment ;
    public final void tree_assigment() throws RecognitionException {
        Token tr=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:69:2: (tr= ( TREE | UTREE ) rich_newick_assignment )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:69:4: tr= ( TREE | UTREE ) rich_newick_assignment
            {
            tr=(Token)input.LT(1);

            if ( input.LA(1)==TREE||input.LA(1)==UTREE ) {
                input.consume();
                state.errorRecovery=false;
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                throw mse;
            }


            pushFollow(FOLLOW_rich_newick_assignment_in_tree_assigment158);
            rich_newick_assignment();

            state._fsp--;


             stack.pushTreeAssignment(tr); 

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
    // $ANTLR end "tree_assigment"



    // $ANTLR start "rich_newick_assignment"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:71:1: rich_newick_assignment : (d= DEFAULT_INDICATOR )? identifier '=' rich_newick_string ;
    public final void rich_newick_assignment() throws RecognitionException {
        Token d=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:72:2: ( (d= DEFAULT_INDICATOR )? identifier '=' rich_newick_string )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:72:4: (d= DEFAULT_INDICATOR )? identifier '=' rich_newick_string
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:72:5: (d= DEFAULT_INDICATOR )?
            int alt6=2;
            int LA6_0 = input.LA(1);

            if ( (LA6_0==DEFAULT_INDICATOR) ) {
                alt6=1;
            }
            switch (alt6) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:72:5: d= DEFAULT_INDICATOR
                    {
                    d=(Token)match(input,DEFAULT_INDICATOR,FOLLOW_DEFAULT_INDICATOR_in_rich_newick_assignment172); 

                    }
                    break;

            }


            pushFollow(FOLLOW_identifier_in_rich_newick_assignment175);
            identifier();

            state._fsp--;


            match(input,19,FOLLOW_19_in_rich_newick_assignment177); 

            pushFollow(FOLLOW_rich_newick_string_in_rich_newick_assignment179);
            rich_newick_string();

            state._fsp--;


             stack.pushRichNewickAssignment(d==null); 

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
    // $ANTLR end "rich_newick_assignment"



    // $ANTLR start "translation"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:74:1: translation : TRANSLATE identifier identifier ( ',' identifier identifier )* ';' ;
    public final void translation() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:75:2: ( TRANSLATE identifier identifier ( ',' identifier identifier )* ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:75:4: TRANSLATE identifier identifier ( ',' identifier identifier )* ';'
            {
            match(input,TRANSLATE,FOLLOW_TRANSLATE_in_translation192); 

            pushFollow(FOLLOW_identifier_in_translation194);
            identifier();

            state._fsp--;


            pushFollow(FOLLOW_identifier_in_translation196);
            identifier();

            state._fsp--;


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:75:36: ( ',' identifier identifier )*
            loop7:
            do {
                int alt7=2;
                int LA7_0 = input.LA(1);

                if ( (LA7_0==17) ) {
                    alt7=1;
                }


                switch (alt7) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:75:37: ',' identifier identifier
            	    {
            	    match(input,17,FOLLOW_17_in_translation199); 

            	    pushFollow(FOLLOW_identifier_in_translation201);
            	    identifier();

            	    state._fsp--;


            	    pushFollow(FOLLOW_identifier_in_translation203);
            	    identifier();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop7;
                }
            } while (true);


            match(input,18,FOLLOW_18_in_translation207); 

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
    // $ANTLR end "translation"



    // $ANTLR start "non_networks_block_body"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:77:1: non_networks_block_body : ~ ( NETWORKS | TREES ) ( . )* ;
    public final void non_networks_block_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:78:2: (~ ( NETWORKS | TREES ) ( . )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:78:4: ~ ( NETWORKS | TREES ) ( . )*
            {
            if ( (input.LA(1) >= BEGIN && input.LA(1) <= NETWORK)||(input.LA(1) >= START && input.LA(1) <= TREE)||(input.LA(1) >= UTREE && input.LA(1) <= 19) ) {
                input.consume();
                state.errorRecovery=false;
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                throw mse;
            }


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:78:24: ( . )*
            loop8:
            do {
                int alt8=2;
                int LA8_0 = input.LA(1);

                if ( (LA8_0==END) ) {
                    alt8=2;
                }
                else if ( ((LA8_0 >= BEGIN && LA8_0 <= DEFAULT_INDICATOR)||(LA8_0 >= ID && LA8_0 <= 19)) ) {
                    alt8=1;
                }


                switch (alt8) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:78:24: .
            	    {
            	    matchAny(input); 

            	    }
            	    break;

            	default :
            	    break loop8;
                }
            } while (true);


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
    // $ANTLR end "non_networks_block_body"



    // $ANTLR start "rich_newick_string"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:81:1: rich_newick_string : str= ( ( . )* ';' ) ;
    public final void rich_newick_string() throws RecognitionException {
        Token str=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:82:2: (str= ( ( . )* ';' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:82:5: str= ( ( . )* ';' )
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:82:9: ( ( . )* ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:82:10: ( . )* ';'
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:82:10: ( . )*
            loop9:
            do {
                int alt9=2;
                int LA9_0 = input.LA(1);

                if ( (LA9_0==18) ) {
                    alt9=2;
                }
                else if ( ((LA9_0 >= BEGIN && LA9_0 <= 17)||LA9_0==19) ) {
                    alt9=1;
                }


                switch (alt9) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:82:10: .
            	    {
            	    matchAny(input); 

            	    }
            	    break;

            	default :
            	    break loop9;
                }
            } while (true);


            match(input,18,FOLLOW_18_in_rich_newick_string246); 

            }


             stack.pushRichNewickString((str!=null?str.getText():null)); 

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
    // $ANTLR end "rich_newick_string"



    // $ANTLR start "identifier"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:84:1: identifier : s= ( TRANSLATE | TREE | UTREE | NETWORK | BEGIN | NETWORKS | TREES | ID ) ;
    public final void identifier() throws RecognitionException {
        Token s=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:85:2: (s= ( TRANSLATE | TREE | UTREE | NETWORK | BEGIN | NETWORKS | TREES | ID ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:85:4: s= ( TRANSLATE | TREE | UTREE | NETWORK | BEGIN | NETWORKS | TREES | ID )
            {
            s=(Token)input.LT(1);

            if ( input.LA(1)==BEGIN||input.LA(1)==ID||(input.LA(1) >= NETWORK && input.LA(1) <= NETWORKS)||(input.LA(1) >= TRANSLATE && input.LA(1) <= UTREE) ) {
                input.consume();
                state.errorRecovery=false;
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                throw mse;
            }


             stack.pushIdentifier((s!=null?s.getText():null)); 

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
    // $ANTLR end "identifier"

    // Delegated rules


 

    public static final BitSet FOLLOW_START_in_blocks36 = new BitSet(new long[]{0x0000000000000012L});
    public static final BitSet FOLLOW_block_in_blocks38 = new BitSet(new long[]{0x0000000000000012L});
    public static final BitSet FOLLOW_BEGIN_in_block49 = new BitSet(new long[]{0x0000000000000400L});
    public static final BitSet FOLLOW_networks_block_body_in_block51 = new BitSet(new long[]{0x0000000000000040L});
    public static final BitSet FOLLOW_END_in_block57 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block62 = new BitSet(new long[]{0x0000000000004000L});
    public static final BitSet FOLLOW_trees_block_body_in_block64 = new BitSet(new long[]{0x0000000000000040L});
    public static final BitSet FOLLOW_END_in_block73 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block78 = new BitSet(new long[]{0x00000000000FBBF0L});
    public static final BitSet FOLLOW_non_networks_block_body_in_block80 = new BitSet(new long[]{0x0000000000000040L});
    public static final BitSet FOLLOW_END_in_block82 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_NETWORKS_in_networks_block_body93 = new BitSet(new long[]{0x0000000000001000L});
    public static final BitSet FOLLOW_translation_in_networks_block_body95 = new BitSet(new long[]{0x0000000000000200L});
    public static final BitSet FOLLOW_NETWORK_in_networks_block_body97 = new BitSet(new long[]{0x000000000000F6B0L});
    public static final BitSet FOLLOW_rich_newick_assignment_in_networks_block_body99 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_TREES_in_trees_block_body110 = new BitSet(new long[]{0x0000000000001000L});
    public static final BitSet FOLLOW_translation_in_trees_block_body112 = new BitSet(new long[]{0x000000000000A002L});
    public static final BitSet FOLLOW_tree_assigment_in_trees_block_body114 = new BitSet(new long[]{0x000000000000A002L});
    public static final BitSet FOLLOW_TREES_in_trees_block_body122 = new BitSet(new long[]{0x000000000000A002L});
    public static final BitSet FOLLOW_tree_assigment_in_trees_block_body136 = new BitSet(new long[]{0x000000000000A002L});
    public static final BitSet FOLLOW_set_in_tree_assigment150 = new BitSet(new long[]{0x000000000000F6B0L});
    public static final BitSet FOLLOW_rich_newick_assignment_in_tree_assigment158 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_DEFAULT_INDICATOR_in_rich_newick_assignment172 = new BitSet(new long[]{0x000000000000F690L});
    public static final BitSet FOLLOW_identifier_in_rich_newick_assignment175 = new BitSet(new long[]{0x0000000000080000L});
    public static final BitSet FOLLOW_19_in_rich_newick_assignment177 = new BitSet(new long[]{0x00000000000FFFF0L});
    public static final BitSet FOLLOW_rich_newick_string_in_rich_newick_assignment179 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_TRANSLATE_in_translation192 = new BitSet(new long[]{0x000000000000F690L});
    public static final BitSet FOLLOW_identifier_in_translation194 = new BitSet(new long[]{0x000000000000F690L});
    public static final BitSet FOLLOW_identifier_in_translation196 = new BitSet(new long[]{0x0000000000060000L});
    public static final BitSet FOLLOW_17_in_translation199 = new BitSet(new long[]{0x000000000000F690L});
    public static final BitSet FOLLOW_identifier_in_translation201 = new BitSet(new long[]{0x000000000000F690L});
    public static final BitSet FOLLOW_identifier_in_translation203 = new BitSet(new long[]{0x0000000000060000L});
    public static final BitSet FOLLOW_18_in_translation207 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_set_in_non_networks_block_body217 = new BitSet(new long[]{0x00000000000FFFF2L});
    public static final BitSet FOLLOW_18_in_rich_newick_string246 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_set_in_identifier261 = new BitSet(new long[]{0x0000000000000002L});

}
// $ANTLR 3.4 D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g 2011-11-03 16:17:38

package edu.rice.cs.bioinfo.library.language.pyson._1_0.parsers.antlr.ast;
import java.util.LinkedList;



import org.antlr.runtime.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;

@SuppressWarnings({"all", "warnings", "unchecked"})
public class PySONParser extends org.antlr.runtime.Parser {
    public static final String[] tokenNames = new String[] {
        "<invalid>", "<EOR>", "<DOWN>", "<UP>", "BEGIN", "DEFAULT_INDICATOR", "ELSE", "END", "ID", "ID_SET", "NESTED_ML_COMMENT", "NETWORK", "NETWORKS", "PHYLONET", "QUOTE", "ROOTAGE_QUALIFIER", "START", "TAXON_SET_LIST", "TRANSLATE", "TREE", "TREES", "UTREE", "WS", "'('", "')'", "','", "':'", "';'", "'<'", "'='", "'>'"
    };

    public static final int EOF=-1;
    public static final int T__23=23;
    public static final int T__24=24;
    public static final int T__25=25;
    public static final int T__26=26;
    public static final int T__27=27;
    public static final int T__28=28;
    public static final int T__29=29;
    public static final int T__30=30;
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
    public static final int TAXON_SET_LIST=17;
    public static final int TRANSLATE=18;
    public static final int TREE=19;
    public static final int TREES=20;
    public static final int UTREE=21;
    public static final int WS=22;

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

    @Override   
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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:1: blocks : START ( block )* ;
    public final void blocks() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:9: ( START ( block )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:11: START ( block )*
            {
            match(input,START,FOLLOW_START_in_blocks52); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:17: ( block )*
            loop1:
            do {
                int alt1=2;
                int LA1_0 = input.LA(1);

                if ( (LA1_0==BEGIN) ) {
                    alt1=1;
                }


                switch (alt1) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:17: block
            	    {
            	    pushFollow(FOLLOW_block_in_blocks54);
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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:91:1: block : ( BEGIN networks_block_body end_semi | BEGIN trees_block_body end_semi | BEGIN phylonet_block_body end_semi | BEGIN skip_body end_semi );
    public final void block() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:91:7: ( BEGIN networks_block_body end_semi | BEGIN trees_block_body end_semi | BEGIN phylonet_block_body end_semi | BEGIN skip_body end_semi )
            int alt2=4;
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
                case PHYLONET:
                    {
                    alt2=3;
                    }
                    break;
                case BEGIN:
                case DEFAULT_INDICATOR:
                case ELSE:
                case END:
                case ID:
                case ID_SET:
                case NESTED_ML_COMMENT:
                case NETWORK:
                case QUOTE:
                case ROOTAGE_QUALIFIER:
                case START:
                case TAXON_SET_LIST:
                case TRANSLATE:
                case TREE:
                case UTREE:
                case WS:
                case 23:
                case 24:
                case 25:
                case 26:
                case 27:
                case 28:
                case 29:
                case 30:
                    {
                    alt2=4;
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
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:91:9: BEGIN networks_block_body end_semi
                    {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block65); 

                    pushFollow(FOLLOW_networks_block_body_in_block67);
                    networks_block_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block69);
                    end_semi();

                    state._fsp--;


                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:92:4: BEGIN trees_block_body end_semi
                    {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block74); 

                    pushFollow(FOLLOW_trees_block_body_in_block76);
                    trees_block_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block81);
                    end_semi();

                    state._fsp--;


                    }
                    break;
                case 3 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:93:4: BEGIN phylonet_block_body end_semi
                    {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block86); 

                    pushFollow(FOLLOW_phylonet_block_body_in_block88);
                    phylonet_block_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block90);
                    end_semi();

                    state._fsp--;


                    }
                    break;
                case 4 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:95:4: BEGIN skip_body end_semi
                    {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block96); 

                    pushFollow(FOLLOW_skip_body_in_block98);
                    skip_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block110);
                    end_semi();

                    state._fsp--;


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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:97:1: networks_block_body : NETWORKS ';' ( NETWORK rich_newick_assignment )* ;
    public final void networks_block_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:98:2: ( NETWORKS ';' ( NETWORK rich_newick_assignment )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:99:3: NETWORKS ';' ( NETWORK rich_newick_assignment )*
            {
            match(input,NETWORKS,FOLLOW_NETWORKS_in_networks_block_body121); 

            match(input,27,FOLLOW_27_in_networks_block_body123); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:99:28: ( NETWORK rich_newick_assignment )*
            loop3:
            do {
                int alt3=2;
                int LA3_0 = input.LA(1);

                if ( (LA3_0==NETWORK) ) {
                    alt3=1;
                }


                switch (alt3) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:99:29: NETWORK rich_newick_assignment
            	    {
            	    match(input,NETWORK,FOLLOW_NETWORK_in_networks_block_body138); 

            	    pushFollow(FOLLOW_rich_newick_assignment_in_networks_block_body140);
            	    rich_newick_assignment();

            	    state._fsp--;


            	     stack.pushNetworkAssignment(); 

            	    }
            	    break;

            	default :
            	    break loop3;
                }
            } while (true);


             stack.pushNetworksBlockBody(false); 

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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:101:1: trees_block_body : TREES ';' ( tree_assigment )* ;
    public final void trees_block_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:103:2: ( TREES ';' ( tree_assigment )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:103:4: TREES ';' ( tree_assigment )*
            {
            match(input,TREES,FOLLOW_TREES_in_trees_block_body160); 

            match(input,27,FOLLOW_27_in_trees_block_body162); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:103:27: ( tree_assigment )*
            loop4:
            do {
                int alt4=2;
                int LA4_0 = input.LA(1);

                if ( (LA4_0==TREE||LA4_0==UTREE) ) {
                    alt4=1;
                }


                switch (alt4) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:103:27: tree_assigment
            	    {
            	    pushFollow(FOLLOW_tree_assigment_in_trees_block_body177);
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



    // $ANTLR start "phylonet_block_body"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:105:1: phylonet_block_body : PHYLONET ';' ( phylonet_command )* ;
    public final void phylonet_block_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:106:2: ( PHYLONET ';' ( phylonet_command )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:106:4: PHYLONET ';' ( phylonet_command )*
            {
            match(input,PHYLONET,FOLLOW_PHYLONET_in_phylonet_block_body189); 

            match(input,27,FOLLOW_27_in_phylonet_block_body191); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:106:17: ( phylonet_command )*
            loop5:
            do {
                int alt5=2;
                int LA5_0 = input.LA(1);

                if ( (LA5_0==BEGIN||(LA5_0 >= ID && LA5_0 <= ID_SET)||(LA5_0 >= NETWORK && LA5_0 <= QUOTE)||(LA5_0 >= TAXON_SET_LIST && LA5_0 <= UTREE)||LA5_0==23||(LA5_0 >= 27 && LA5_0 <= 28)) ) {
                    alt5=1;
                }


                switch (alt5) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:106:17: phylonet_command
            	    {
            	    pushFollow(FOLLOW_phylonet_command_in_phylonet_block_body193);
            	    phylonet_command();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop5;
                }
            } while (true);


             stack.pushPhylonetBlockBody(); 

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
    // $ANTLR end "phylonet_block_body"



    // $ANTLR start "phylonet_command"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:114:1: phylonet_command : ( phylonet_command_part )* ';' ;
    public final void phylonet_command() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:115:2: ( ( phylonet_command_part )* ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:115:5: ( phylonet_command_part )* ';'
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:115:5: ( phylonet_command_part )*
            loop6:
            do {
                int alt6=2;
                int LA6_0 = input.LA(1);

                if ( (LA6_0==BEGIN||(LA6_0 >= ID && LA6_0 <= ID_SET)||(LA6_0 >= NETWORK && LA6_0 <= QUOTE)||(LA6_0 >= TAXON_SET_LIST && LA6_0 <= UTREE)||LA6_0==23||LA6_0==28) ) {
                    alt6=1;
                }


                switch (alt6) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:115:5: phylonet_command_part
            	    {
            	    pushFollow(FOLLOW_phylonet_command_part_in_phylonet_command214);
            	    phylonet_command_part();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop6;
                }
            } while (true);


            match(input,27,FOLLOW_27_in_phylonet_command217); 

             stack.pushPhylonetCommand(); 

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
    // $ANTLR end "phylonet_command"



    // $ANTLR start "phylonet_command_part"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:117:1: phylonet_command_part : ( identifier | ident_list |p= QUOTE |p= TAXON_SET_LIST |p= ID_SET | taxa_map );
    public final void phylonet_command_part() throws RecognitionException {
        Token p=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:118:2: ( identifier | ident_list |p= QUOTE |p= TAXON_SET_LIST |p= ID_SET | taxa_map )
            int alt7=6;
            switch ( input.LA(1) ) {
            case BEGIN:
            case ID:
            case NETWORK:
            case NETWORKS:
            case PHYLONET:
            case TRANSLATE:
            case TREE:
            case TREES:
            case UTREE:
                {
                alt7=1;
                }
                break;
            case 23:
                {
                alt7=2;
                }
                break;
            case QUOTE:
                {
                alt7=3;
                }
                break;
            case TAXON_SET_LIST:
                {
                alt7=4;
                }
                break;
            case ID_SET:
                {
                alt7=5;
                }
                break;
            case 28:
                {
                alt7=6;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 7, 0, input);

                throw nvae;

            }

            switch (alt7) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:118:4: identifier
                    {
                    pushFollow(FOLLOW_identifier_in_phylonet_command_part229);
                    identifier();

                    state._fsp--;


                     stack.pushPhylonetCommandPartIdent(); 

                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:119:4: ident_list
                    {
                    pushFollow(FOLLOW_ident_list_in_phylonet_command_part248);
                    ident_list();

                    state._fsp--;


                     stack.pushPhylonetCommandPartIdentList(); 

                    }
                    break;
                case 3 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:120:4: p= QUOTE
                    {
                    p=(Token)match(input,QUOTE,FOLLOW_QUOTE_in_phylonet_command_part264); 

                     stack.pushPhylonetCommandPartQuote(p); 

                    }
                    break;
                case 4 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:121:10: p= TAXON_SET_LIST
                    {
                    p=(Token)match(input,TAXON_SET_LIST,FOLLOW_TAXON_SET_LIST_in_phylonet_command_part294); 

                     stack.pushPhylonetCommandPartSetList(p); 

                    }
                    break;
                case 5 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:122:4: p= ID_SET
                    {
                    p=(Token)match(input,ID_SET,FOLLOW_ID_SET_in_phylonet_command_part309); 

                     stack.pushPhylonetCommandPartIdSet(p); 

                    }
                    break;
                case 6 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:123:4: taxa_map
                    {
                    pushFollow(FOLLOW_taxa_map_in_phylonet_command_part323);
                    taxa_map();

                    state._fsp--;


                     stack.pushPhylonetCommandPartTaxaMap(); 

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
    // $ANTLR end "phylonet_command_part"



    // $ANTLR start "ident_list"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:125:1: ident_list : s= '(' identifier ( ',' identifier )* ')' ;
    public final void ident_list() throws RecognitionException {
        Token s=null;

         int numIdentsInList = 1; 
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:126:2: (s= '(' identifier ( ',' identifier )* ')' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:126:4: s= '(' identifier ( ',' identifier )* ')'
            {
            s=(Token)match(input,23,FOLLOW_23_in_ident_list355); 

            pushFollow(FOLLOW_identifier_in_ident_list357);
            identifier();

            state._fsp--;


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:126:21: ( ',' identifier )*
            loop8:
            do {
                int alt8=2;
                int LA8_0 = input.LA(1);

                if ( (LA8_0==25) ) {
                    alt8=1;
                }


                switch (alt8) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:126:22: ',' identifier
            	    {
            	    match(input,25,FOLLOW_25_in_ident_list360); 

            	    pushFollow(FOLLOW_identifier_in_ident_list362);
            	    identifier();

            	    state._fsp--;


            	     numIdentsInList++; 

            	    }
            	    break;

            	default :
            	    break loop8;
                }
            } while (true);


            match(input,24,FOLLOW_24_in_ident_list369); 

             stack.pushIdentList(numIdentsInList, s); 

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
    // $ANTLR end "ident_list"



    // $ANTLR start "tree_assigment"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:128:1: tree_assigment : tr= ( TREE | UTREE ) rich_newick_assignment ;
    public final void tree_assigment() throws RecognitionException {
        Token tr=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:129:2: (tr= ( TREE | UTREE ) rich_newick_assignment )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:129:4: tr= ( TREE | UTREE ) rich_newick_assignment
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


            pushFollow(FOLLOW_rich_newick_assignment_in_tree_assigment391);
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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:131:1: rich_newick_assignment : (d= DEFAULT_INDICATOR )? identifier '=' rich_newick_string ;
    public final void rich_newick_assignment() throws RecognitionException {
        Token d=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:132:2: ( (d= DEFAULT_INDICATOR )? identifier '=' rich_newick_string )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:132:4: (d= DEFAULT_INDICATOR )? identifier '=' rich_newick_string
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:132:5: (d= DEFAULT_INDICATOR )?
            int alt9=2;
            int LA9_0 = input.LA(1);

            if ( (LA9_0==DEFAULT_INDICATOR) ) {
                alt9=1;
            }
            switch (alt9) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:132:5: d= DEFAULT_INDICATOR
                    {
                    d=(Token)match(input,DEFAULT_INDICATOR,FOLLOW_DEFAULT_INDICATOR_in_rich_newick_assignment405); 

                    }
                    break;

            }


            pushFollow(FOLLOW_identifier_in_rich_newick_assignment408);
            identifier();

            state._fsp--;


            match(input,29,FOLLOW_29_in_rich_newick_assignment410); 

            pushFollow(FOLLOW_rich_newick_string_in_rich_newick_assignment412);
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



    // $ANTLR start "skip_body"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:135:1: skip_body : ~ ( TREES | NETWORKS | PHYLONET ) ';' (~ END )* ;
    public final void skip_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:136:2: (~ ( TREES | NETWORKS | PHYLONET ) ';' (~ END )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:136:4: ~ ( TREES | NETWORKS | PHYLONET ) ';' (~ END )*
            {
            if ( (input.LA(1) >= BEGIN && input.LA(1) <= NETWORK)||(input.LA(1) >= QUOTE && input.LA(1) <= TREE)||(input.LA(1) >= UTREE && input.LA(1) <= 30) ) {
                input.consume();
                state.errorRecovery=false;
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                throw mse;
            }


            match(input,27,FOLLOW_27_in_skip_body435); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:136:35: (~ END )*
            loop10:
            do {
                int alt10=2;
                int LA10_0 = input.LA(1);

                if ( ((LA10_0 >= BEGIN && LA10_0 <= ELSE)||(LA10_0 >= ID && LA10_0 <= 30)) ) {
                    alt10=1;
                }


                switch (alt10) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
            	    {
            	    if ( (input.LA(1) >= BEGIN && input.LA(1) <= ELSE)||(input.LA(1) >= ID && input.LA(1) <= 30) ) {
            	        input.consume();
            	        state.errorRecovery=false;
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        throw mse;
            	    }


            	    }
            	    break;

            	default :
            	    break loop10;
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
    // $ANTLR end "skip_body"



    // $ANTLR start "rich_newick_string"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:139:1: rich_newick_string : (str=~ ( ';' ) )* ';' ;
    public final void rich_newick_string() throws RecognitionException {
        Token str=null;

         StringBuffer accum = new StringBuffer(); int line = -1; int col = -1; 
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:140:2: ( (str=~ ( ';' ) )* ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:140:5: (str=~ ( ';' ) )* ';'
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:140:5: (str=~ ( ';' ) )*
            loop11:
            do {
                int alt11=2;
                int LA11_0 = input.LA(1);

                if ( ((LA11_0 >= BEGIN && LA11_0 <= 26)||(LA11_0 >= 28 && LA11_0 <= 30)) ) {
                    alt11=1;
                }


                switch (alt11) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:140:6: str=~ ( ';' )
            	    {
            	    str=(Token)input.LT(1);

            	    if ( (input.LA(1) >= BEGIN && input.LA(1) <= 26)||(input.LA(1) >= 28 && input.LA(1) <= 30) ) {
            	        input.consume();
            	        state.errorRecovery=false;
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        throw mse;
            	    }


            	    accum.append((str!=null?str.getText():null)); 
            	    	                      line = line == -1 ? (str!=null?str.getLine():0)    : line; 
            	    	                      col  = col  == -1 ? (str!=null?str.getCharPositionInLine():0)     : col ; 

            	    }
            	    break;

            	default :
            	    break loop11;
                }
            } while (true);


            match(input,27,FOLLOW_27_in_rich_newick_string481); 

             stack.pushRichNewickString(accum.toString() + ';', line, col); 

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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:145:1: identifier : s= ( TRANSLATE | TREE | UTREE | NETWORK | BEGIN | NETWORKS | TREES | PHYLONET | ID ) ;
    public final void identifier() throws RecognitionException {
        Token s=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:146:2: (s= ( TRANSLATE | TREE | UTREE | NETWORK | BEGIN | NETWORKS | TREES | PHYLONET | ID ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:146:4: s= ( TRANSLATE | TREE | UTREE | NETWORK | BEGIN | NETWORKS | TREES | PHYLONET | ID )
            {
            s=(Token)input.LT(1);

            if ( input.LA(1)==BEGIN||input.LA(1)==ID||(input.LA(1) >= NETWORK && input.LA(1) <= PHYLONET)||(input.LA(1) >= TRANSLATE && input.LA(1) <= UTREE) ) {
                input.consume();
                state.errorRecovery=false;
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                throw mse;
            }


             stack.pushIdentifier(s); 

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



    // $ANTLR start "end_semi"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:148:1: end_semi : END ';' ;
    public final void end_semi() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:148:9: ( END ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:148:11: END ';'
            {
            match(input,END,FOLLOW_END_in_end_semi540); 

            match(input,27,FOLLOW_27_in_end_semi542); 

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
    // $ANTLR end "end_semi"



    // $ANTLR start "taxa_map"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:150:1: taxa_map : s= '<' ( taxa_map_entry ';' )+ '>' ;
    public final void taxa_map() throws RecognitionException {
        Token s=null;

         int numKeys = 0; 
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:150:36: (s= '<' ( taxa_map_entry ';' )+ '>' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:152:2: s= '<' ( taxa_map_entry ';' )+ '>'
            {
            s=(Token)match(input,28,FOLLOW_28_in_taxa_map560); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:152:7: ( taxa_map_entry ';' )+
            int cnt12=0;
            loop12:
            do {
                int alt12=2;
                int LA12_0 = input.LA(1);

                if ( (LA12_0==BEGIN||LA12_0==ID||(LA12_0 >= NETWORK && LA12_0 <= PHYLONET)||(LA12_0 >= TRANSLATE && LA12_0 <= UTREE)) ) {
                    alt12=1;
                }


                switch (alt12) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:152:9: taxa_map_entry ';'
            	    {
            	    pushFollow(FOLLOW_taxa_map_entry_in_taxa_map563);
            	    taxa_map_entry();

            	    state._fsp--;


            	    match(input,27,FOLLOW_27_in_taxa_map565); 

            	    numKeys++;

            	    }
            	    break;

            	default :
            	    if ( cnt12 >= 1 ) break loop12;
                        EarlyExitException eee =
                            new EarlyExitException(12, input);
                        throw eee;
                }
                cnt12++;
            } while (true);


            match(input,30,FOLLOW_30_in_taxa_map570); 

             stack.pushTaxaMap(numKeys, s); 

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
    // $ANTLR end "taxa_map"



    // $ANTLR start "taxa_map_entry"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:1: taxa_map_entry : identifier ':' identifier ( ',' identifier )* ;
    public final void taxa_map_entry() throws RecognitionException {
         int numValues = 1; 
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:155:2: ( identifier ':' identifier ( ',' identifier )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:155:4: identifier ':' identifier ( ',' identifier )*
            {
            pushFollow(FOLLOW_identifier_in_taxa_map_entry586);
            identifier();

            state._fsp--;


            match(input,26,FOLLOW_26_in_taxa_map_entry588); 

            pushFollow(FOLLOW_identifier_in_taxa_map_entry590);
            identifier();

            state._fsp--;


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:155:30: ( ',' identifier )*
            loop13:
            do {
                int alt13=2;
                int LA13_0 = input.LA(1);

                if ( (LA13_0==25) ) {
                    alt13=1;
                }


                switch (alt13) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:155:31: ',' identifier
            	    {
            	    match(input,25,FOLLOW_25_in_taxa_map_entry593); 

            	    pushFollow(FOLLOW_identifier_in_taxa_map_entry595);
            	    identifier();

            	    state._fsp--;


            	     numValues++; 

            	    }
            	    break;

            	default :
            	    break loop13;
                }
            } while (true);


             stack.pushTaxaMapEntry(numValues); 

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
    // $ANTLR end "taxa_map_entry"

    // Delegated rules


 

    public static final BitSet FOLLOW_START_in_blocks52 = new BitSet(new long[]{0x0000000000000012L});
    public static final BitSet FOLLOW_block_in_blocks54 = new BitSet(new long[]{0x0000000000000012L});
    public static final BitSet FOLLOW_BEGIN_in_block65 = new BitSet(new long[]{0x0000000000001000L});
    public static final BitSet FOLLOW_networks_block_body_in_block67 = new BitSet(new long[]{0x0000000000000080L});
    public static final BitSet FOLLOW_end_semi_in_block69 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block74 = new BitSet(new long[]{0x0000000000100000L});
    public static final BitSet FOLLOW_trees_block_body_in_block76 = new BitSet(new long[]{0x0000000000000080L});
    public static final BitSet FOLLOW_end_semi_in_block81 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block86 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_phylonet_block_body_in_block88 = new BitSet(new long[]{0x0000000000000080L});
    public static final BitSet FOLLOW_end_semi_in_block90 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block96 = new BitSet(new long[]{0x000000007FEFCFF0L});
    public static final BitSet FOLLOW_skip_body_in_block98 = new BitSet(new long[]{0x0000000000000080L});
    public static final BitSet FOLLOW_end_semi_in_block110 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_NETWORKS_in_networks_block_body121 = new BitSet(new long[]{0x0000000008000000L});
    public static final BitSet FOLLOW_27_in_networks_block_body123 = new BitSet(new long[]{0x0000000000000802L});
    public static final BitSet FOLLOW_NETWORK_in_networks_block_body138 = new BitSet(new long[]{0x00000000003C3930L});
    public static final BitSet FOLLOW_rich_newick_assignment_in_networks_block_body140 = new BitSet(new long[]{0x0000000000000802L});
    public static final BitSet FOLLOW_TREES_in_trees_block_body160 = new BitSet(new long[]{0x0000000008000000L});
    public static final BitSet FOLLOW_27_in_trees_block_body162 = new BitSet(new long[]{0x0000000000280002L});
    public static final BitSet FOLLOW_tree_assigment_in_trees_block_body177 = new BitSet(new long[]{0x0000000000280002L});
    public static final BitSet FOLLOW_PHYLONET_in_phylonet_block_body189 = new BitSet(new long[]{0x0000000008000000L});
    public static final BitSet FOLLOW_27_in_phylonet_block_body191 = new BitSet(new long[]{0x0000000018BE7B12L});
    public static final BitSet FOLLOW_phylonet_command_in_phylonet_block_body193 = new BitSet(new long[]{0x0000000018BE7B12L});
    public static final BitSet FOLLOW_phylonet_command_part_in_phylonet_command214 = new BitSet(new long[]{0x0000000018BE7B10L});
    public static final BitSet FOLLOW_27_in_phylonet_command217 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_identifier_in_phylonet_command_part229 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_ident_list_in_phylonet_command_part248 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_QUOTE_in_phylonet_command_part264 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_TAXON_SET_LIST_in_phylonet_command_part294 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_ID_SET_in_phylonet_command_part309 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_taxa_map_in_phylonet_command_part323 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_23_in_ident_list355 = new BitSet(new long[]{0x00000000003C3910L});
    public static final BitSet FOLLOW_identifier_in_ident_list357 = new BitSet(new long[]{0x0000000003000000L});
    public static final BitSet FOLLOW_25_in_ident_list360 = new BitSet(new long[]{0x00000000003C3910L});
    public static final BitSet FOLLOW_identifier_in_ident_list362 = new BitSet(new long[]{0x0000000003000000L});
    public static final BitSet FOLLOW_24_in_ident_list369 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_set_in_tree_assigment383 = new BitSet(new long[]{0x00000000003C3930L});
    public static final BitSet FOLLOW_rich_newick_assignment_in_tree_assigment391 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_DEFAULT_INDICATOR_in_rich_newick_assignment405 = new BitSet(new long[]{0x00000000003C3910L});
    public static final BitSet FOLLOW_identifier_in_rich_newick_assignment408 = new BitSet(new long[]{0x0000000020000000L});
    public static final BitSet FOLLOW_29_in_rich_newick_assignment410 = new BitSet(new long[]{0x000000007FFFFFF0L});
    public static final BitSet FOLLOW_rich_newick_string_in_rich_newick_assignment412 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_set_in_skip_body426 = new BitSet(new long[]{0x0000000008000000L});
    public static final BitSet FOLLOW_27_in_skip_body435 = new BitSet(new long[]{0x000000007FFFFF72L});
    public static final BitSet FOLLOW_set_in_rich_newick_string460 = new BitSet(new long[]{0x000000007FFFFFF0L});
    public static final BitSet FOLLOW_27_in_rich_newick_string481 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_set_in_identifier495 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_END_in_end_semi540 = new BitSet(new long[]{0x0000000008000000L});
    public static final BitSet FOLLOW_27_in_end_semi542 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_28_in_taxa_map560 = new BitSet(new long[]{0x00000000003C3910L});
    public static final BitSet FOLLOW_taxa_map_entry_in_taxa_map563 = new BitSet(new long[]{0x0000000008000000L});
    public static final BitSet FOLLOW_27_in_taxa_map565 = new BitSet(new long[]{0x00000000403C3910L});
    public static final BitSet FOLLOW_30_in_taxa_map570 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_identifier_in_taxa_map_entry586 = new BitSet(new long[]{0x0000000004000000L});
    public static final BitSet FOLLOW_26_in_taxa_map_entry588 = new BitSet(new long[]{0x00000000003C3910L});
    public static final BitSet FOLLOW_identifier_in_taxa_map_entry590 = new BitSet(new long[]{0x0000000002000002L});
    public static final BitSet FOLLOW_25_in_taxa_map_entry593 = new BitSet(new long[]{0x00000000003C3910L});
    public static final BitSet FOLLOW_identifier_in_taxa_map_entry595 = new BitSet(new long[]{0x0000000002000002L});

}
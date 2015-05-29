// $ANTLR 3.4 D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g 2015-05-29 14:53:08

package edu.rice.cs.bioinfo.library.language.pyson._1_0.parsers.antlr.ast;
import java.util.LinkedList;



import org.antlr.runtime.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;

@SuppressWarnings({"all", "warnings", "unchecked"})
public class PySONParser extends org.antlr.runtime.Parser {
    public static final String[] tokenNames = new String[] {
        "<invalid>", "<EOR>", "<DOWN>", "<UP>", "BEGIN", "DATA", "DATATYPE", "DEFAULT_INDICATOR", "DIMENSIONS", "ELSE", "END", "FORMAT", "GAP", "ID", "ID_SET", "MATRIX", "MISSING", "MORPHDATA", "NCHAR", "NESTED_ML_COMMENT", "NETWORK", "NETWORKS", "NTAX", "PHYLONET", "QUOTE", "RN_LS_NONCOMMENT", "START", "SYMBOLS", "TAXON_SET_LIST", "TRANSLATE", "TREE", "TREES", "UTREE", "WS", "'('", "')'", "','", "':'", "';'", "'<'", "'='", "'>'"
    };

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
    public String getGrammarFileName() { return "D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g"; }



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
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:1: blocks : START ( block )* EOF ;
    public final void blocks() throws RecognitionException {
        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:9: ( START ( block )* EOF )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:11: START ( block )* EOF
            {
            match(input,START,FOLLOW_START_in_blocks52); 

            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:17: ( block )*
            loop1:
            do {
                int alt1=2;
                int LA1_0 = input.LA(1);

                if ( (LA1_0==BEGIN) ) {
                    alt1=1;
                }


                switch (alt1) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:17: block
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

            match(input,EOF,FOLLOW_EOF_in_blocks59); 

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
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:91:1: block : ( BEGIN networks_block_body end_semi | BEGIN trees_block_body end_semi | BEGIN phylonet_block_body end_semi | BEGIN data_block_body end_semi | BEGIN morph_block_body end_semi | BEGIN skip_body end_semi );
    public final void block() throws RecognitionException {
        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:91:7: ( BEGIN networks_block_body end_semi | BEGIN trees_block_body end_semi | BEGIN phylonet_block_body end_semi | BEGIN data_block_body end_semi | BEGIN morph_block_body end_semi | BEGIN skip_body end_semi )
            int alt2=6;
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
                case DATA:
                    {
                    alt2=4;
                    }
                    break;
                case MORPHDATA:
                    {
                    alt2=5;
                    }
                    break;
                case BEGIN:
                case DATATYPE:
                case DEFAULT_INDICATOR:
                case DIMENSIONS:
                case ELSE:
                case END:
                case FORMAT:
                case GAP:
                case ID:
                case ID_SET:
                case MATRIX:
                case MISSING:
                case NCHAR:
                case NESTED_ML_COMMENT:
                case NETWORK:
                case NTAX:
                case QUOTE:
                case RN_LS_NONCOMMENT:
                case START:
                case SYMBOLS:
                case TAXON_SET_LIST:
                case TRANSLATE:
                case TREE:
                case UTREE:
                case WS:
                case 34:
                case 35:
                case 36:
                case 37:
                case 38:
                case 39:
                case 40:
                case 41:
                    {
                    alt2=6;
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
                    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:91:9: BEGIN networks_block_body end_semi
                    {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block67); 

                    pushFollow(FOLLOW_networks_block_body_in_block69);
                    networks_block_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block71);
                    end_semi();

                    state._fsp--;


                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:92:4: BEGIN trees_block_body end_semi
                    {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block76); 

                    pushFollow(FOLLOW_trees_block_body_in_block78);
                    trees_block_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block83);
                    end_semi();

                    state._fsp--;


                    }
                    break;
                case 3 :
                    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:93:4: BEGIN phylonet_block_body end_semi
                    {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block88); 

                    pushFollow(FOLLOW_phylonet_block_body_in_block90);
                    phylonet_block_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block92);
                    end_semi();

                    state._fsp--;


                    }
                    break;
                case 4 :
                    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:94:4: BEGIN data_block_body end_semi
                    {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block97); 

                    pushFollow(FOLLOW_data_block_body_in_block99);
                    data_block_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block103);
                    end_semi();

                    state._fsp--;


                    }
                    break;
                case 5 :
                    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:95:4: BEGIN morph_block_body end_semi
                    {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block108); 

                    pushFollow(FOLLOW_morph_block_body_in_block110);
                    morph_block_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block114);
                    end_semi();

                    state._fsp--;


                    }
                    break;
                case 6 :
                    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:96:4: BEGIN skip_body end_semi
                    {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block119); 

                    pushFollow(FOLLOW_skip_body_in_block121);
                    skip_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block133);
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
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:98:1: networks_block_body : NETWORKS ';' ( NETWORK rich_newick_assignment )* ;
    public final void networks_block_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:99:2: ( NETWORKS ';' ( NETWORK rich_newick_assignment )* )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:100:3: NETWORKS ';' ( NETWORK rich_newick_assignment )*
            {
            match(input,NETWORKS,FOLLOW_NETWORKS_in_networks_block_body144); 

            match(input,38,FOLLOW_38_in_networks_block_body146); 

            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:100:28: ( NETWORK rich_newick_assignment )*
            loop3:
            do {
                int alt3=2;
                int LA3_0 = input.LA(1);

                if ( (LA3_0==NETWORK) ) {
                    alt3=1;
                }


                switch (alt3) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:100:29: NETWORK rich_newick_assignment
            	    {
            	    match(input,NETWORK,FOLLOW_NETWORK_in_networks_block_body161); 

            	    pushFollow(FOLLOW_rich_newick_assignment_in_networks_block_body163);
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
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:102:1: trees_block_body : TREES ';' ( tree_assigment )* ;
    public final void trees_block_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:104:2: ( TREES ';' ( tree_assigment )* )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:104:4: TREES ';' ( tree_assigment )*
            {
            match(input,TREES,FOLLOW_TREES_in_trees_block_body183); 

            match(input,38,FOLLOW_38_in_trees_block_body185); 

            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:104:27: ( tree_assigment )*
            loop4:
            do {
                int alt4=2;
                int LA4_0 = input.LA(1);

                if ( (LA4_0==TREE||LA4_0==UTREE) ) {
                    alt4=1;
                }


                switch (alt4) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:104:27: tree_assigment
            	    {
            	    pushFollow(FOLLOW_tree_assigment_in_trees_block_body200);
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
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:106:1: phylonet_block_body : PHYLONET ';' ( phylonet_command )* ;
    public final void phylonet_block_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:107:2: ( PHYLONET ';' ( phylonet_command )* )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:107:4: PHYLONET ';' ( phylonet_command )*
            {
            match(input,PHYLONET,FOLLOW_PHYLONET_in_phylonet_block_body212); 

            match(input,38,FOLLOW_38_in_phylonet_block_body214); 

            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:107:17: ( phylonet_command )*
            loop5:
            do {
                int alt5=2;
                int LA5_0 = input.LA(1);

                if ( ((LA5_0 >= BEGIN && LA5_0 <= DATA)||LA5_0==DIMENSIONS||(LA5_0 >= FORMAT && LA5_0 <= MISSING)||LA5_0==NCHAR||(LA5_0 >= NETWORK && LA5_0 <= QUOTE)||(LA5_0 >= SYMBOLS && LA5_0 <= UTREE)||LA5_0==34||(LA5_0 >= 38 && LA5_0 <= 39)) ) {
                    alt5=1;
                }


                switch (alt5) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:107:17: phylonet_command
            	    {
            	    pushFollow(FOLLOW_phylonet_command_in_phylonet_block_body216);
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



    // $ANTLR start "data_block_body"
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:109:1: data_block_body : DATA ';' ( DIMENSIONS NTAX '=' identifier NCHAR '=' identifier ';' )? FORMAT DATATYPE '=' identifier SYMBOLS '=' QUOTE MISSING '=' ID GAP '=' ID ';' MATRIX ( identifier identifier )* ';' ;
    public final void data_block_body() throws RecognitionException {
         int numPairs = 0; 
        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:110:2: ( DATA ';' ( DIMENSIONS NTAX '=' identifier NCHAR '=' identifier ';' )? FORMAT DATATYPE '=' identifier SYMBOLS '=' QUOTE MISSING '=' ID GAP '=' ID ';' MATRIX ( identifier identifier )* ';' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:110:4: DATA ';' ( DIMENSIONS NTAX '=' identifier NCHAR '=' identifier ';' )? FORMAT DATATYPE '=' identifier SYMBOLS '=' QUOTE MISSING '=' ID GAP '=' ID ';' MATRIX ( identifier identifier )* ';'
            {
            match(input,DATA,FOLLOW_DATA_in_data_block_body233); 

            match(input,38,FOLLOW_38_in_data_block_body235); 

            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:111:3: ( DIMENSIONS NTAX '=' identifier NCHAR '=' identifier ';' )?
            int alt6=2;
            int LA6_0 = input.LA(1);

            if ( (LA6_0==DIMENSIONS) ) {
                alt6=1;
            }
            switch (alt6) {
                case 1 :
                    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:111:4: DIMENSIONS NTAX '=' identifier NCHAR '=' identifier ';'
                    {
                    match(input,DIMENSIONS,FOLLOW_DIMENSIONS_in_data_block_body240); 

                    match(input,NTAX,FOLLOW_NTAX_in_data_block_body242); 

                    match(input,40,FOLLOW_40_in_data_block_body244); 

                    pushFollow(FOLLOW_identifier_in_data_block_body246);
                    identifier();

                    state._fsp--;


                    match(input,NCHAR,FOLLOW_NCHAR_in_data_block_body248); 

                    match(input,40,FOLLOW_40_in_data_block_body250); 

                    pushFollow(FOLLOW_identifier_in_data_block_body252);
                    identifier();

                    state._fsp--;


                    match(input,38,FOLLOW_38_in_data_block_body254); 

                    }
                    break;

            }


            match(input,FORMAT,FOLLOW_FORMAT_in_data_block_body260); 

            match(input,DATATYPE,FOLLOW_DATATYPE_in_data_block_body262); 

            match(input,40,FOLLOW_40_in_data_block_body263); 

            pushFollow(FOLLOW_identifier_in_data_block_body264);
            identifier();

            state._fsp--;


            match(input,SYMBOLS,FOLLOW_SYMBOLS_in_data_block_body266); 

            match(input,40,FOLLOW_40_in_data_block_body267); 

            match(input,QUOTE,FOLLOW_QUOTE_in_data_block_body268); 

            match(input,MISSING,FOLLOW_MISSING_in_data_block_body270); 

            match(input,40,FOLLOW_40_in_data_block_body271); 

            match(input,ID,FOLLOW_ID_in_data_block_body272); 

            match(input,GAP,FOLLOW_GAP_in_data_block_body274); 

            match(input,40,FOLLOW_40_in_data_block_body275); 

            match(input,ID,FOLLOW_ID_in_data_block_body276); 

            match(input,38,FOLLOW_38_in_data_block_body278); 

            match(input,MATRIX,FOLLOW_MATRIX_in_data_block_body282); 

            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:114:3: ( identifier identifier )*
            loop7:
            do {
                int alt7=2;
                int LA7_0 = input.LA(1);

                if ( ((LA7_0 >= BEGIN && LA7_0 <= DATA)||LA7_0==DIMENSIONS||(LA7_0 >= FORMAT && LA7_0 <= ID)||(LA7_0 >= MATRIX && LA7_0 <= MISSING)||LA7_0==NCHAR||(LA7_0 >= NETWORK && LA7_0 <= PHYLONET)||LA7_0==SYMBOLS||(LA7_0 >= TRANSLATE && LA7_0 <= UTREE)) ) {
                    alt7=1;
                }


                switch (alt7) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:114:4: identifier identifier
            	    {
            	    pushFollow(FOLLOW_identifier_in_data_block_body287);
            	    identifier();

            	    state._fsp--;


            	    pushFollow(FOLLOW_identifier_in_data_block_body289);
            	    identifier();

            	    state._fsp--;


            	    numPairs++;

            	    }
            	    break;

            	default :
            	    break loop7;
                }
            } while (true);


            stack.pushDataBlockBody(numPairs); 

            match(input,38,FOLLOW_38_in_data_block_body298); 

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
    // $ANTLR end "data_block_body"



    // $ANTLR start "morph_block_body"
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:116:1: morph_block_body : MORPHDATA ';' ( DIMENSIONS NTAX '=' identifier )? FORMAT SYMBOLS '=' s= QUOTE MISSING '=' m= ID ';' MATRIX ( identifier )* ';' ;
    public final void morph_block_body() throws RecognitionException {
        Token s=null;
        Token m=null;

         int numIdents = 0; 
        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:117:2: ( MORPHDATA ';' ( DIMENSIONS NTAX '=' identifier )? FORMAT SYMBOLS '=' s= QUOTE MISSING '=' m= ID ';' MATRIX ( identifier )* ';' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:117:4: MORPHDATA ';' ( DIMENSIONS NTAX '=' identifier )? FORMAT SYMBOLS '=' s= QUOTE MISSING '=' m= ID ';' MATRIX ( identifier )* ';'
            {
            match(input,MORPHDATA,FOLLOW_MORPHDATA_in_morph_block_body312); 

            match(input,38,FOLLOW_38_in_morph_block_body314); 

            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:118:3: ( DIMENSIONS NTAX '=' identifier )?
            int alt8=2;
            int LA8_0 = input.LA(1);

            if ( (LA8_0==DIMENSIONS) ) {
                alt8=1;
            }
            switch (alt8) {
                case 1 :
                    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:118:4: DIMENSIONS NTAX '=' identifier
                    {
                    match(input,DIMENSIONS,FOLLOW_DIMENSIONS_in_morph_block_body319); 

                    match(input,NTAX,FOLLOW_NTAX_in_morph_block_body321); 

                    match(input,40,FOLLOW_40_in_morph_block_body323); 

                    pushFollow(FOLLOW_identifier_in_morph_block_body325);
                    identifier();

                    state._fsp--;


                    }
                    break;

            }


            match(input,FORMAT,FOLLOW_FORMAT_in_morph_block_body331); 

            match(input,SYMBOLS,FOLLOW_SYMBOLS_in_morph_block_body333); 

            match(input,40,FOLLOW_40_in_morph_block_body334); 

            s=(Token)match(input,QUOTE,FOLLOW_QUOTE_in_morph_block_body337); 

            match(input,MISSING,FOLLOW_MISSING_in_morph_block_body339); 

            match(input,40,FOLLOW_40_in_morph_block_body340); 

            m=(Token)match(input,ID,FOLLOW_ID_in_morph_block_body343); 

            match(input,38,FOLLOW_38_in_morph_block_body344); 

            match(input,MATRIX,FOLLOW_MATRIX_in_morph_block_body348); 

            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:121:3: ( identifier )*
            loop9:
            do {
                int alt9=2;
                int LA9_0 = input.LA(1);

                if ( ((LA9_0 >= BEGIN && LA9_0 <= DATA)||LA9_0==DIMENSIONS||(LA9_0 >= FORMAT && LA9_0 <= ID)||(LA9_0 >= MATRIX && LA9_0 <= MISSING)||LA9_0==NCHAR||(LA9_0 >= NETWORK && LA9_0 <= PHYLONET)||LA9_0==SYMBOLS||(LA9_0 >= TRANSLATE && LA9_0 <= UTREE)) ) {
                    alt9=1;
                }


                switch (alt9) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:121:4: identifier
            	    {
            	    pushFollow(FOLLOW_identifier_in_morph_block_body353);
            	    identifier();

            	    state._fsp--;


            	    numIdents++;

            	    }
            	    break;

            	default :
            	    break loop9;
                }
            } while (true);


            stack.pushMorphDataBlockBody(numIdents, s.getText(), m.getText());

            match(input,38,FOLLOW_38_in_morph_block_body362); 

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
    // $ANTLR end "morph_block_body"



    // $ANTLR start "phylonet_command"
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:126:1: phylonet_command : ( identifier '=' )? ( phylonet_command_part )* ';' ;
    public final void phylonet_command() throws RecognitionException {
         boolean assignment = false; 
        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:127:2: ( ( identifier '=' )? ( phylonet_command_part )* ';' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:127:5: ( identifier '=' )? ( phylonet_command_part )* ';'
            {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:127:5: ( identifier '=' )?
            int alt10=2;
            int LA10_0 = input.LA(1);

            if ( ((LA10_0 >= BEGIN && LA10_0 <= DATA)||LA10_0==DIMENSIONS||(LA10_0 >= FORMAT && LA10_0 <= ID)||(LA10_0 >= MATRIX && LA10_0 <= MISSING)||LA10_0==NCHAR||(LA10_0 >= NETWORK && LA10_0 <= PHYLONET)||LA10_0==SYMBOLS||(LA10_0 >= TRANSLATE && LA10_0 <= UTREE)) ) {
                int LA10_1 = input.LA(2);

                if ( (LA10_1==40) ) {
                    alt10=1;
                }
            }
            switch (alt10) {
                case 1 :
                    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:127:6: identifier '='
                    {
                    pushFollow(FOLLOW_identifier_in_phylonet_command393);
                    identifier();

                    state._fsp--;


                    match(input,40,FOLLOW_40_in_phylonet_command395); 

                    assignment = true;

                    }
                    break;

            }


            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:127:45: ( phylonet_command_part )*
            loop11:
            do {
                int alt11=2;
                int LA11_0 = input.LA(1);

                if ( ((LA11_0 >= BEGIN && LA11_0 <= DATA)||LA11_0==DIMENSIONS||(LA11_0 >= FORMAT && LA11_0 <= MISSING)||LA11_0==NCHAR||(LA11_0 >= NETWORK && LA11_0 <= QUOTE)||(LA11_0 >= SYMBOLS && LA11_0 <= UTREE)||LA11_0==34||LA11_0==39) ) {
                    alt11=1;
                }


                switch (alt11) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:127:45: phylonet_command_part
            	    {
            	    pushFollow(FOLLOW_phylonet_command_part_in_phylonet_command402);
            	    phylonet_command_part();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop11;
                }
            } while (true);


            match(input,38,FOLLOW_38_in_phylonet_command405); 

             stack.pushPhylonetCommand(assignment); 

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
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:129:1: phylonet_command_part : ( identifier | ident_list |p= QUOTE |p= TAXON_SET_LIST |p= ID_SET | taxa_map );
    public final void phylonet_command_part() throws RecognitionException {
        Token p=null;

        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:130:2: ( identifier | ident_list |p= QUOTE |p= TAXON_SET_LIST |p= ID_SET | taxa_map )
            int alt12=6;
            switch ( input.LA(1) ) {
            case BEGIN:
            case DATA:
            case DIMENSIONS:
            case FORMAT:
            case GAP:
            case ID:
            case MATRIX:
            case MISSING:
            case NCHAR:
            case NETWORK:
            case NETWORKS:
            case NTAX:
            case PHYLONET:
            case SYMBOLS:
            case TRANSLATE:
            case TREE:
            case TREES:
            case UTREE:
                {
                alt12=1;
                }
                break;
            case 34:
                {
                alt12=2;
                }
                break;
            case QUOTE:
                {
                alt12=3;
                }
                break;
            case TAXON_SET_LIST:
                {
                alt12=4;
                }
                break;
            case ID_SET:
                {
                alt12=5;
                }
                break;
            case 39:
                {
                alt12=6;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 12, 0, input);

                throw nvae;

            }

            switch (alt12) {
                case 1 :
                    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:130:4: identifier
                    {
                    pushFollow(FOLLOW_identifier_in_phylonet_command_part417);
                    identifier();

                    state._fsp--;


                     stack.pushPhylonetCommandPartIdent(); 

                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:131:4: ident_list
                    {
                    pushFollow(FOLLOW_ident_list_in_phylonet_command_part436);
                    ident_list();

                    state._fsp--;


                     stack.pushPhylonetCommandPartIdentList(); 

                    }
                    break;
                case 3 :
                    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:132:4: p= QUOTE
                    {
                    p=(Token)match(input,QUOTE,FOLLOW_QUOTE_in_phylonet_command_part452); 

                     stack.pushPhylonetCommandPartQuote(p); 

                    }
                    break;
                case 4 :
                    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:133:10: p= TAXON_SET_LIST
                    {
                    p=(Token)match(input,TAXON_SET_LIST,FOLLOW_TAXON_SET_LIST_in_phylonet_command_part482); 

                     stack.pushPhylonetCommandPartSetList(p); 

                    }
                    break;
                case 5 :
                    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:134:4: p= ID_SET
                    {
                    p=(Token)match(input,ID_SET,FOLLOW_ID_SET_in_phylonet_command_part497); 

                     stack.pushPhylonetCommandPartIdSet(p); 

                    }
                    break;
                case 6 :
                    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:135:4: taxa_map
                    {
                    pushFollow(FOLLOW_taxa_map_in_phylonet_command_part511);
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
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:137:1: ident_list : s= '(' identifier ( ',' identifier )* ')' ;
    public final void ident_list() throws RecognitionException {
        Token s=null;

         int numIdentsInList = 1; 
        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:138:2: (s= '(' identifier ( ',' identifier )* ')' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:138:4: s= '(' identifier ( ',' identifier )* ')'
            {
            s=(Token)match(input,34,FOLLOW_34_in_ident_list543); 

            pushFollow(FOLLOW_identifier_in_ident_list545);
            identifier();

            state._fsp--;


            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:138:21: ( ',' identifier )*
            loop13:
            do {
                int alt13=2;
                int LA13_0 = input.LA(1);

                if ( (LA13_0==36) ) {
                    alt13=1;
                }


                switch (alt13) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:138:22: ',' identifier
            	    {
            	    match(input,36,FOLLOW_36_in_ident_list548); 

            	    pushFollow(FOLLOW_identifier_in_ident_list550);
            	    identifier();

            	    state._fsp--;


            	     numIdentsInList++; 

            	    }
            	    break;

            	default :
            	    break loop13;
                }
            } while (true);


            match(input,35,FOLLOW_35_in_ident_list557); 

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
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:140:1: tree_assigment : tr= ( TREE | UTREE ) rich_newick_assignment ;
    public final void tree_assigment() throws RecognitionException {
        Token tr=null;

        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:141:2: (tr= ( TREE | UTREE ) rich_newick_assignment )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:141:4: tr= ( TREE | UTREE ) rich_newick_assignment
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


            pushFollow(FOLLOW_rich_newick_assignment_in_tree_assigment579);
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
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:143:1: rich_newick_assignment : (d= DEFAULT_INDICATOR )? identifier '=' rich_newick_string ;
    public final void rich_newick_assignment() throws RecognitionException {
        Token d=null;

        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:144:2: ( (d= DEFAULT_INDICATOR )? identifier '=' rich_newick_string )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:144:4: (d= DEFAULT_INDICATOR )? identifier '=' rich_newick_string
            {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:144:5: (d= DEFAULT_INDICATOR )?
            int alt14=2;
            int LA14_0 = input.LA(1);

            if ( (LA14_0==DEFAULT_INDICATOR) ) {
                alt14=1;
            }
            switch (alt14) {
                case 1 :
                    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:144:5: d= DEFAULT_INDICATOR
                    {
                    d=(Token)match(input,DEFAULT_INDICATOR,FOLLOW_DEFAULT_INDICATOR_in_rich_newick_assignment593); 

                    }
                    break;

            }


            pushFollow(FOLLOW_identifier_in_rich_newick_assignment596);
            identifier();

            state._fsp--;


            match(input,40,FOLLOW_40_in_rich_newick_assignment598); 

            pushFollow(FOLLOW_rich_newick_string_in_rich_newick_assignment600);
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
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:147:1: skip_body : ~ ( TREES | NETWORKS | PHYLONET | DATA | MORPHDATA ) ';' (~ END )* ;
    public final void skip_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:148:2: (~ ( TREES | NETWORKS | PHYLONET | DATA | MORPHDATA ) ';' (~ END )* )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:148:4: ~ ( TREES | NETWORKS | PHYLONET | DATA | MORPHDATA ) ';' (~ END )*
            {
            if ( input.LA(1)==BEGIN||(input.LA(1) >= DATATYPE && input.LA(1) <= MISSING)||(input.LA(1) >= NCHAR && input.LA(1) <= NETWORK)||input.LA(1)==NTAX||(input.LA(1) >= QUOTE && input.LA(1) <= TREE)||(input.LA(1) >= UTREE && input.LA(1) <= 41) ) {
                input.consume();
                state.errorRecovery=false;
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                throw mse;
            }


            match(input,38,FOLLOW_38_in_skip_body627); 

            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:148:50: (~ END )*
            loop15:
            do {
                int alt15=2;
                int LA15_0 = input.LA(1);

                if ( ((LA15_0 >= BEGIN && LA15_0 <= ELSE)||(LA15_0 >= FORMAT && LA15_0 <= 41)) ) {
                    alt15=1;
                }


                switch (alt15) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
            	    {
            	    if ( (input.LA(1) >= BEGIN && input.LA(1) <= ELSE)||(input.LA(1) >= FORMAT && input.LA(1) <= 41) ) {
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
            	    break loop15;
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
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:151:1: rich_newick_string : (str=~ ( ';' ) )* ';' ;
    public final void rich_newick_string() throws RecognitionException {
        Token str=null;

         StringBuffer accum = new StringBuffer(); int line = -1; int col = -1; 
        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:152:2: ( (str=~ ( ';' ) )* ';' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:152:5: (str=~ ( ';' ) )* ';'
            {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:152:5: (str=~ ( ';' ) )*
            loop16:
            do {
                int alt16=2;
                int LA16_0 = input.LA(1);

                if ( ((LA16_0 >= BEGIN && LA16_0 <= 37)||(LA16_0 >= 39 && LA16_0 <= 41)) ) {
                    alt16=1;
                }


                switch (alt16) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:152:6: str=~ ( ';' )
            	    {
            	    str=(Token)input.LT(1);

            	    if ( (input.LA(1) >= BEGIN && input.LA(1) <= 37)||(input.LA(1) >= 39 && input.LA(1) <= 41) ) {
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
            	    break loop16;
                }
            } while (true);


            match(input,38,FOLLOW_38_in_rich_newick_string673); 

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
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:157:1: identifier : s= ( TRANSLATE | DATA | DIMENSIONS | SYMBOLS | MISSING | GAP | MATRIX | TREE | UTREE | FORMAT | NETWORK | NCHAR | BEGIN | NETWORKS | NTAX | TREES | PHYLONET | ID ) ;
    public final void identifier() throws RecognitionException {
        Token s=null;

        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:158:2: (s= ( TRANSLATE | DATA | DIMENSIONS | SYMBOLS | MISSING | GAP | MATRIX | TREE | UTREE | FORMAT | NETWORK | NCHAR | BEGIN | NETWORKS | NTAX | TREES | PHYLONET | ID ) )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:158:4: s= ( TRANSLATE | DATA | DIMENSIONS | SYMBOLS | MISSING | GAP | MATRIX | TREE | UTREE | FORMAT | NETWORK | NCHAR | BEGIN | NETWORKS | NTAX | TREES | PHYLONET | ID )
            {
            s=(Token)input.LT(1);

            if ( (input.LA(1) >= BEGIN && input.LA(1) <= DATA)||input.LA(1)==DIMENSIONS||(input.LA(1) >= FORMAT && input.LA(1) <= ID)||(input.LA(1) >= MATRIX && input.LA(1) <= MISSING)||input.LA(1)==NCHAR||(input.LA(1) >= NETWORK && input.LA(1) <= PHYLONET)||input.LA(1)==SYMBOLS||(input.LA(1) >= TRANSLATE && input.LA(1) <= UTREE) ) {
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
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:160:1: end_semi : END ';' ;
    public final void end_semi() throws RecognitionException {
        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:160:9: ( END ';' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:160:11: END ';'
            {
            match(input,END,FOLLOW_END_in_end_semi768); 

            match(input,38,FOLLOW_38_in_end_semi770); 

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
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:162:1: taxa_map : s= '<' taxa_map_entry ( ';' taxa_map_entry )+ '>' ;
    public final void taxa_map() throws RecognitionException {
        Token s=null;

         int numKeys = 1; 
        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:162:36: (s= '<' taxa_map_entry ( ';' taxa_map_entry )+ '>' )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:164:2: s= '<' taxa_map_entry ( ';' taxa_map_entry )+ '>'
            {
            s=(Token)match(input,39,FOLLOW_39_in_taxa_map788); 

            pushFollow(FOLLOW_taxa_map_entry_in_taxa_map789);
            taxa_map_entry();

            state._fsp--;


            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:164:22: ( ';' taxa_map_entry )+
            int cnt17=0;
            loop17:
            do {
                int alt17=2;
                int LA17_0 = input.LA(1);

                if ( (LA17_0==38) ) {
                    alt17=1;
                }


                switch (alt17) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:164:23: ';' taxa_map_entry
            	    {
            	    match(input,38,FOLLOW_38_in_taxa_map792); 

            	    pushFollow(FOLLOW_taxa_map_entry_in_taxa_map794);
            	    taxa_map_entry();

            	    state._fsp--;


            	    numKeys++;

            	    }
            	    break;

            	default :
            	    if ( cnt17 >= 1 ) break loop17;
                        EarlyExitException eee =
                            new EarlyExitException(17, input);
                        throw eee;
                }
                cnt17++;
            } while (true);


            match(input,41,FOLLOW_41_in_taxa_map800); 

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
    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:166:1: taxa_map_entry : identifier ':' identifier ( ',' identifier )* ;
    public final void taxa_map_entry() throws RecognitionException {
         int numValues = 1; 
        try {
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:167:2: ( identifier ':' identifier ( ',' identifier )* )
            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:167:4: identifier ':' identifier ( ',' identifier )*
            {
            pushFollow(FOLLOW_identifier_in_taxa_map_entry816);
            identifier();

            state._fsp--;


            match(input,37,FOLLOW_37_in_taxa_map_entry818); 

            pushFollow(FOLLOW_identifier_in_taxa_map_entry820);
            identifier();

            state._fsp--;


            // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:167:30: ( ',' identifier )*
            loop18:
            do {
                int alt18=2;
                int LA18_0 = input.LA(1);

                if ( (LA18_0==36) ) {
                    alt18=1;
                }


                switch (alt18) {
            	case 1 :
            	    // D:\\WorkDev\\Bioinfo\\Code\\Antlr\\Unstable\\PySON\\PySON.g:167:31: ',' identifier
            	    {
            	    match(input,36,FOLLOW_36_in_taxa_map_entry823); 

            	    pushFollow(FOLLOW_identifier_in_taxa_map_entry825);
            	    identifier();

            	    state._fsp--;


            	     numValues++; 

            	    }
            	    break;

            	default :
            	    break loop18;
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


 

    public static final BitSet FOLLOW_START_in_blocks52 = new BitSet(new long[]{0x0000000000000010L});
    public static final BitSet FOLLOW_block_in_blocks54 = new BitSet(new long[]{0x0000000000000010L});
    public static final BitSet FOLLOW_EOF_in_blocks59 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block67 = new BitSet(new long[]{0x0000000000200000L});
    public static final BitSet FOLLOW_networks_block_body_in_block69 = new BitSet(new long[]{0x0000000000000400L});
    public static final BitSet FOLLOW_end_semi_in_block71 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block76 = new BitSet(new long[]{0x0000000080000000L});
    public static final BitSet FOLLOW_trees_block_body_in_block78 = new BitSet(new long[]{0x0000000000000400L});
    public static final BitSet FOLLOW_end_semi_in_block83 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block88 = new BitSet(new long[]{0x0000000000800000L});
    public static final BitSet FOLLOW_phylonet_block_body_in_block90 = new BitSet(new long[]{0x0000000000000400L});
    public static final BitSet FOLLOW_end_semi_in_block92 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block97 = new BitSet(new long[]{0x0000000000000020L});
    public static final BitSet FOLLOW_data_block_body_in_block99 = new BitSet(new long[]{0x0000000000000400L});
    public static final BitSet FOLLOW_end_semi_in_block103 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block108 = new BitSet(new long[]{0x0000000000020000L});
    public static final BitSet FOLLOW_morph_block_body_in_block110 = new BitSet(new long[]{0x0000000000000400L});
    public static final BitSet FOLLOW_end_semi_in_block114 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block119 = new BitSet(new long[]{0x000003FF7F5DFFD0L});
    public static final BitSet FOLLOW_skip_body_in_block121 = new BitSet(new long[]{0x0000000000000400L});
    public static final BitSet FOLLOW_end_semi_in_block133 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_NETWORKS_in_networks_block_body144 = new BitSet(new long[]{0x0000004000000000L});
    public static final BitSet FOLLOW_38_in_networks_block_body146 = new BitSet(new long[]{0x0000000000100002L});
    public static final BitSet FOLLOW_NETWORK_in_networks_block_body161 = new BitSet(new long[]{0x00000001E8F5B9B0L});
    public static final BitSet FOLLOW_rich_newick_assignment_in_networks_block_body163 = new BitSet(new long[]{0x0000000000100002L});
    public static final BitSet FOLLOW_TREES_in_trees_block_body183 = new BitSet(new long[]{0x0000004000000000L});
    public static final BitSet FOLLOW_38_in_trees_block_body185 = new BitSet(new long[]{0x0000000140000002L});
    public static final BitSet FOLLOW_tree_assigment_in_trees_block_body200 = new BitSet(new long[]{0x0000000140000002L});
    public static final BitSet FOLLOW_PHYLONET_in_phylonet_block_body212 = new BitSet(new long[]{0x0000004000000000L});
    public static final BitSet FOLLOW_38_in_phylonet_block_body214 = new BitSet(new long[]{0x000000C5F9F5F932L});
    public static final BitSet FOLLOW_phylonet_command_in_phylonet_block_body216 = new BitSet(new long[]{0x000000C5F9F5F932L});
    public static final BitSet FOLLOW_DATA_in_data_block_body233 = new BitSet(new long[]{0x0000004000000000L});
    public static final BitSet FOLLOW_38_in_data_block_body235 = new BitSet(new long[]{0x0000000000000900L});
    public static final BitSet FOLLOW_DIMENSIONS_in_data_block_body240 = new BitSet(new long[]{0x0000000000400000L});
    public static final BitSet FOLLOW_NTAX_in_data_block_body242 = new BitSet(new long[]{0x0000010000000000L});
    public static final BitSet FOLLOW_40_in_data_block_body244 = new BitSet(new long[]{0x00000001E8F5B930L});
    public static final BitSet FOLLOW_identifier_in_data_block_body246 = new BitSet(new long[]{0x0000000000040000L});
    public static final BitSet FOLLOW_NCHAR_in_data_block_body248 = new BitSet(new long[]{0x0000010000000000L});
    public static final BitSet FOLLOW_40_in_data_block_body250 = new BitSet(new long[]{0x00000001E8F5B930L});
    public static final BitSet FOLLOW_identifier_in_data_block_body252 = new BitSet(new long[]{0x0000004000000000L});
    public static final BitSet FOLLOW_38_in_data_block_body254 = new BitSet(new long[]{0x0000000000000800L});
    public static final BitSet FOLLOW_FORMAT_in_data_block_body260 = new BitSet(new long[]{0x0000000000000040L});
    public static final BitSet FOLLOW_DATATYPE_in_data_block_body262 = new BitSet(new long[]{0x0000010000000000L});
    public static final BitSet FOLLOW_40_in_data_block_body263 = new BitSet(new long[]{0x00000001E8F5B930L});
    public static final BitSet FOLLOW_identifier_in_data_block_body264 = new BitSet(new long[]{0x0000000008000000L});
    public static final BitSet FOLLOW_SYMBOLS_in_data_block_body266 = new BitSet(new long[]{0x0000010000000000L});
    public static final BitSet FOLLOW_40_in_data_block_body267 = new BitSet(new long[]{0x0000000001000000L});
    public static final BitSet FOLLOW_QUOTE_in_data_block_body268 = new BitSet(new long[]{0x0000000000010000L});
    public static final BitSet FOLLOW_MISSING_in_data_block_body270 = new BitSet(new long[]{0x0000010000000000L});
    public static final BitSet FOLLOW_40_in_data_block_body271 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_ID_in_data_block_body272 = new BitSet(new long[]{0x0000000000001000L});
    public static final BitSet FOLLOW_GAP_in_data_block_body274 = new BitSet(new long[]{0x0000010000000000L});
    public static final BitSet FOLLOW_40_in_data_block_body275 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_ID_in_data_block_body276 = new BitSet(new long[]{0x0000004000000000L});
    public static final BitSet FOLLOW_38_in_data_block_body278 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_MATRIX_in_data_block_body282 = new BitSet(new long[]{0x00000041E8F5B930L});
    public static final BitSet FOLLOW_identifier_in_data_block_body287 = new BitSet(new long[]{0x00000001E8F5B930L});
    public static final BitSet FOLLOW_identifier_in_data_block_body289 = new BitSet(new long[]{0x00000041E8F5B930L});
    public static final BitSet FOLLOW_38_in_data_block_body298 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_MORPHDATA_in_morph_block_body312 = new BitSet(new long[]{0x0000004000000000L});
    public static final BitSet FOLLOW_38_in_morph_block_body314 = new BitSet(new long[]{0x0000000000000900L});
    public static final BitSet FOLLOW_DIMENSIONS_in_morph_block_body319 = new BitSet(new long[]{0x0000000000400000L});
    public static final BitSet FOLLOW_NTAX_in_morph_block_body321 = new BitSet(new long[]{0x0000010000000000L});
    public static final BitSet FOLLOW_40_in_morph_block_body323 = new BitSet(new long[]{0x00000001E8F5B930L});
    public static final BitSet FOLLOW_identifier_in_morph_block_body325 = new BitSet(new long[]{0x0000000000000800L});
    public static final BitSet FOLLOW_FORMAT_in_morph_block_body331 = new BitSet(new long[]{0x0000000008000000L});
    public static final BitSet FOLLOW_SYMBOLS_in_morph_block_body333 = new BitSet(new long[]{0x0000010000000000L});
    public static final BitSet FOLLOW_40_in_morph_block_body334 = new BitSet(new long[]{0x0000000001000000L});
    public static final BitSet FOLLOW_QUOTE_in_morph_block_body337 = new BitSet(new long[]{0x0000000000010000L});
    public static final BitSet FOLLOW_MISSING_in_morph_block_body339 = new BitSet(new long[]{0x0000010000000000L});
    public static final BitSet FOLLOW_40_in_morph_block_body340 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_ID_in_morph_block_body343 = new BitSet(new long[]{0x0000004000000000L});
    public static final BitSet FOLLOW_38_in_morph_block_body344 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_MATRIX_in_morph_block_body348 = new BitSet(new long[]{0x00000041E8F5B930L});
    public static final BitSet FOLLOW_identifier_in_morph_block_body353 = new BitSet(new long[]{0x00000041E8F5B930L});
    public static final BitSet FOLLOW_38_in_morph_block_body362 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_identifier_in_phylonet_command393 = new BitSet(new long[]{0x0000010000000000L});
    public static final BitSet FOLLOW_40_in_phylonet_command395 = new BitSet(new long[]{0x000000C5F9F5F930L});
    public static final BitSet FOLLOW_phylonet_command_part_in_phylonet_command402 = new BitSet(new long[]{0x000000C5F9F5F930L});
    public static final BitSet FOLLOW_38_in_phylonet_command405 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_identifier_in_phylonet_command_part417 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_ident_list_in_phylonet_command_part436 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_QUOTE_in_phylonet_command_part452 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_TAXON_SET_LIST_in_phylonet_command_part482 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_ID_SET_in_phylonet_command_part497 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_taxa_map_in_phylonet_command_part511 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_34_in_ident_list543 = new BitSet(new long[]{0x00000001E8F5B930L});
    public static final BitSet FOLLOW_identifier_in_ident_list545 = new BitSet(new long[]{0x0000001800000000L});
    public static final BitSet FOLLOW_36_in_ident_list548 = new BitSet(new long[]{0x00000001E8F5B930L});
    public static final BitSet FOLLOW_identifier_in_ident_list550 = new BitSet(new long[]{0x0000001800000000L});
    public static final BitSet FOLLOW_35_in_ident_list557 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_set_in_tree_assigment571 = new BitSet(new long[]{0x00000001E8F5B9B0L});
    public static final BitSet FOLLOW_rich_newick_assignment_in_tree_assigment579 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_DEFAULT_INDICATOR_in_rich_newick_assignment593 = new BitSet(new long[]{0x00000001E8F5B930L});
    public static final BitSet FOLLOW_identifier_in_rich_newick_assignment596 = new BitSet(new long[]{0x0000010000000000L});
    public static final BitSet FOLLOW_40_in_rich_newick_assignment598 = new BitSet(new long[]{0x000003FFFFFFFFF0L});
    public static final BitSet FOLLOW_rich_newick_string_in_rich_newick_assignment600 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_set_in_skip_body614 = new BitSet(new long[]{0x0000004000000000L});
    public static final BitSet FOLLOW_38_in_skip_body627 = new BitSet(new long[]{0x000003FFFFFFFBF2L});
    public static final BitSet FOLLOW_set_in_rich_newick_string652 = new BitSet(new long[]{0x000003FFFFFFFFF0L});
    public static final BitSet FOLLOW_38_in_rich_newick_string673 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_set_in_identifier687 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_END_in_end_semi768 = new BitSet(new long[]{0x0000004000000000L});
    public static final BitSet FOLLOW_38_in_end_semi770 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_39_in_taxa_map788 = new BitSet(new long[]{0x00000001E8F5B930L});
    public static final BitSet FOLLOW_taxa_map_entry_in_taxa_map789 = new BitSet(new long[]{0x0000004000000000L});
    public static final BitSet FOLLOW_38_in_taxa_map792 = new BitSet(new long[]{0x00000001E8F5B930L});
    public static final BitSet FOLLOW_taxa_map_entry_in_taxa_map794 = new BitSet(new long[]{0x0000024000000000L});
    public static final BitSet FOLLOW_41_in_taxa_map800 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_identifier_in_taxa_map_entry816 = new BitSet(new long[]{0x0000002000000000L});
    public static final BitSet FOLLOW_37_in_taxa_map_entry818 = new BitSet(new long[]{0x00000001E8F5B930L});
    public static final BitSet FOLLOW_identifier_in_taxa_map_entry820 = new BitSet(new long[]{0x0000001000000002L});
    public static final BitSet FOLLOW_36_in_taxa_map_entry823 = new BitSet(new long[]{0x00000001E8F5B930L});
    public static final BitSet FOLLOW_identifier_in_taxa_map_entry825 = new BitSet(new long[]{0x0000001000000002L});

}
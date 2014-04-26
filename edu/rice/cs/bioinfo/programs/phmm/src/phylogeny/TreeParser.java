/**
 * This file is part of PhyloNet-HMM.
 *
 * Copyright © 2013-2014 Kevin Liu, Jingxuan Dai, Kathy Truong, 
 * Ying Song, Michael H. Kohn, and Luay Nakhleh. <http://bioinfo.cs.rice.edu/>
 * 
 * PhyloNet-HMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * PhyloNet-HMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package phylogeny;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.StreamTokenizer;
import java.util.ArrayList;
import java.util.EmptyStackException;
import java.util.Stack;

import javax.swing.JProgressBar;


/**
 * @author James
 *
 * Parses the newick portion of a file
 * For nexus files, additional node-number mapping is needed to rename files
 * Identification of a file as either newick or nexus determines contents
 *
 * */
public class TreeParser
{
    /** Nexus file identifier.  We look for this as the first token to identify a tree file as Nexus, or other. */
    private static final String nexusFileID = "#NEXUS";
    /** Begin tag. */
    private static final String beginTag = "begin";
    /** End tag. */
    private static final String endTag = "end";
    //  trees section
    /** Tree section. */
    private static final String treeSectionTag = "trees";
    /** Tree ID. */
    private static final String treeID = "tree";
    /** Tree ID (same or similar to {@link #treeID}?). */
    private static final String utreeID = "utree"; // two different tree IDs?

    /** Line (and tree information) termination. */
    private static final char lineTerminator = ';';
    /** Equality sign. */
    private static final char equals = '=';
    /** Nexus comment open. */
    private static final char commentOpen = '[';
    /** Nexus comment close. */
    private static final char commentClose = ']';

    /**
     * True: show debug output.  False: suppress printing.
     */
    private StreamTokenizer tokenizer;
    /**
     * Root node of the tree being parsed.  Must be initialized outside the tokenizer.
     */
    private Node rootNode;


    /**
     * Parses names of trees in nexus file.
     * @param fileName Name of nexus file.
     * @return List of all tree names found in nexus file
     */
    @SuppressWarnings({ "rawtypes", "unchecked" })
    public ArrayList<EvoTree> nexusFileTreeNames(String fileName)
    {
        ArrayList returnList = null;
        ArrayList<EvoTree> thetrees = new ArrayList<EvoTree>();
        //BufferedReader r;
        try
        {
            //r = new BufferedReader(new FileReader(fileName));
            StreamTokenizer st = tokenizer;
            st.wordChars('#', '#');
            st.nextToken();
            returnList = new ArrayList();
            while (st.ttype != StreamTokenizer.TT_EOF)
            {
                if (st.ttype == StreamTokenizer.TT_WORD)
                {
                    if (st.sval.equalsIgnoreCase(beginTag))
                    {
                        st.nextToken();
                        if (st.ttype == StreamTokenizer.TT_WORD &&
                                st.sval.equalsIgnoreCase(treeSectionTag))
                        {
                            // found a tree section, huzzah
                            boolean endOfTreeList = false;
                            st.nextToken();
                            while (st.ttype != StreamTokenizer.TT_EOF && !endOfTreeList)
                            {
                                // expect either a tree/utree id or the end tag
                                if (st.ttype == StreamTokenizer.TT_WORD)
                                {
                                    if (st.sval.equalsIgnoreCase(endTag))
                                        endOfTreeList = true;
                                    else if (st.sval.equalsIgnoreCase(treeID) ||
                                            st.sval.equalsIgnoreCase(utreeID))
                                    {
                                        // found the start of a tree
                                        st.nextToken();
                                        if (st.ttype == StreamTokenizer.TT_WORD)
                                        {
                                            returnList.add(st.sval); // found a tree name
                                        }
                                        thetrees.add(tokenize(1, st.sval, null));
//                                        while (st.nextToken() != StreamTokenizer.TT_EOF &&
//                                                st.ttype != ';'); // find the end of the tree
                                    }
                                }
                                else st.nextToken(); // eat a non-word while looking for first tree word

//                                    System.out.println("Not a word while looking for a tree start tag: " + st.ttype);
                            }
                        }
                        // not a tree section, find the end tag or the next start tag
                        else while (st.nextToken() != StreamTokenizer.TT_EOF &&
                                    st.ttype != StreamTokenizer.TT_WORD ||
                                    (!st.sval.equalsIgnoreCase(beginTag) &&
                                    !st.sval.equalsIgnoreCase(endTag)));
                    }
                    else
                        st.nextToken();
                }
                else
                    st.nextToken();
            }
            //r.close();
        }
        catch (FileNotFoundException e)
        {
            System.err.println("Could not find file to identify: " + fileName);
        }
        catch (IOException e)
        {
            System.out.println("Couldn't identify file: " + fileName);
        }
        return thetrees;
    }


    /**
     * Initializes parsing of a tree by creating a tokenizer and setting default
     * properties (such as spacing, quoting characters).
     * {@link #tokenize(long, String, JProgressBar)} is required to start the parsing.
     * @param b Buffered reader that could start in the middle of a nexus file or
     * the start of a newick file (basically the beginning of a newick tree, is run
     * for each tree in a nexus file)
     */
    public TreeParser(BufferedReader b)
    {
        tokenizer = new StreamTokenizer(b);
        tokenizer.eolIsSignificant(false);
        tokenizer.quoteChar('"');
//        tokenizer.quoteChar('\''); // TODO: check quote layering, quoted quotes
        tokenizer.wordChars('\'', '\''); // quote problem, turn this into a prime symbol?
        // 32 = space
        tokenizer.wordChars('!', '!'); // 33
        // 34 = "
        tokenizer.wordChars('#', '&'); // 35-38
        // 39-41 = '() newick
        tokenizer.wordChars('*', '+'); // 42-43
        // 44 = , newick
        tokenizer.wordChars('-', '/'); // 45-47
        // 48-59 = [0-9]:;
        tokenizer.wordChars('<', '<'); // 60
        // 61 = = nexus
        tokenizer.wordChars('>', '@'); // 62-64
        // 65-90 = [A-Z]
//        tokenizer.wordChars('[', '['); // 91 [ nexus comment character, treat as char
        // 92 = \ (esc, support esc'd spaces)
//      93 = ] nexus comment character
        tokenizer.wordChars('^', '`'); // 93-96
        // 97-122 = [a-z]
        tokenizer.wordChars('{', '~'); // 123-126
        // 127 = del
    }

    /**
     * Adds node at the top of the stack to the tree.  TreeNode is already created based
     * on Newick properties.
     * @param name Name of the node.
     * @param nodeStack Stack of nodes that haven't been added to the tree yet.  Nodes are popped when
     * they have names and all children are processed.
     * @return Newly added treeNode linked into the tree.
     */
    private Node popAndName(String name, @SuppressWarnings("rawtypes") Stack nodeStack)
    {
        Node topNode = (Node)nodeStack.pop();
        if (name == null)
        {
            //
        }
        else
        {
            topNode.setTaxa(name);
        }
        try
        {
            Node parent = (Node) nodeStack.peek();
            topNode.setParent(parent);
            parent.addChild(topNode);
        }
        catch (EmptyStackException e)
        {
            if (topNode != rootNode)
                System.out.println("Parser error on node " + topNode);
        }
        //topNode.setExtremeLeaves(); // sets leftmost and rightmost leaf, non-recursive
        //topNode.setNumberLeaves(); // sets number of leaves, non-recursive
        //topNode.linkNodesInPreorder();
        //topNode.linkNodesInPostorder();
        return topNode;
    }

    /**
     * Newick tokenizer: converts a string (tree as a string) into a tree object.
     * The stream tokenizer should be initialized before calling this function.
     * @param fileLength Length of the file, for progress bar movements.
     * For nexus files, this would be the relative position of the next semicolon = the size of the tree in bytes.
     * @param streamName Name of the tree or file that is being loaded.  Nexus files have names ("tree <name> = ((...));", newick trees are named by file name.
     * @param progressBar Reference to a progress bar widgit, embedded perhaps in place of the new canvas for this tree.  If this is null, create a new progress bar here.
     * @return Tree parsed from the stream.
     */
    @SuppressWarnings("unchecked")
    public EvoTree tokenize(long fileLength, String streamName,
            JProgressBar progressBar)
    {
        final char openBracket = '(', closeBracket = ')', childSeparator = ',',
            treeTerminator = lineTerminator, quote = '\'', doubleQuote = '"', infoSeparator = ':';
        int progress = 0;
        rootNode = new Node();
        EvoTree t = new EvoTree();
        t.setRoot(rootNode);
        t.setName(streamName);
        @SuppressWarnings("rawtypes")
        Stack nodeStack = new Stack();
        nodeStack.push(rootNode);
        int thisToken;
        Node lastNamed = null;
        boolean EOT = false;
        boolean nameNext = true;
        int percentage = 0;
    try {
            while (EOT == false &&
                    (thisToken = tokenizer.nextToken()) != StreamTokenizer.TT_EOF)
            {
            switch (thisToken)
            {
//            	case quote:
                case doubleQuote:
                case StreamTokenizer.TT_WORD:
                    if (!nameNext)
                        System.err.println("Error: didn't expect this name here: " + tokenizer.sval);
                    lastNamed = popAndName(tokenizer.sval, nodeStack);
                    progress += tokenizer.sval.length();
                    nameNext = false;
                    break;
                case StreamTokenizer.TT_NUMBER:
                    if (nameNext)
                        lastNamed = popAndName(tokenizer.sval, nodeStack);
                    else
                    {
                        if (lastNamed != null)
                            lastNamed.setTbranch(tokenizer.nval);
                        else
                            System.err.println("Error: can't set value " + tokenizer.nval + " to a null node");
                        lastNamed = null;
                    }
                    progress += (new Double(tokenizer.nval).toString()).length();
                    nameNext = false;
                    break;
                case infoSeparator:
                    if (nameNext)
                        lastNamed = popAndName(null, nodeStack);
                    progress += 1;
                    nameNext = false;
                    break;
                case treeTerminator:
                case StreamTokenizer.TT_EOF:
                    if (nameNext)
                        lastNamed = popAndName(null, nodeStack);
                    EOT = true;
                    progress += 1;
                    nameNext = false;
                    break;
                case openBracket:
                    nodeStack.push(new Node());
                    progress += 1;
                    nameNext = true;
                    break;
                case closeBracket:
                    if (nameNext)
                        lastNamed = popAndName(null, nodeStack);
                    progress += 1;
                    nameNext = true;
                    break;
                case childSeparator:
                    if (nameNext)
                        lastNamed = popAndName(null, nodeStack);
                    nodeStack.push(new Node());
                    progress += 1;
                    nameNext = true;
                    break;
                default:
                    //debugOutput("default " + (char)thisToken);
                    break;
            }
        }
        }
        catch (IOException e) {
        }
        if (!nodeStack.isEmpty())
            System.err.println("Node stack still has " + nodeStack.size() + " things");
        //t.postProcess();
        return t;
    }


}

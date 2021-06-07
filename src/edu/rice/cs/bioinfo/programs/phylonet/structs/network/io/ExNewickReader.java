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

package edu.rice.cs.bioinfo.programs.phylonet.structs.network.io;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;

import java.io.*;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;

/**
 * Provides methods for reading a network in the extended Newick (eNewick) format.
 *
 * @param <T> 	the object type that parsed networks will contain
 */
public class ExNewickReader<T> {
	// Constructors

	/**
	 * Instantiates a new network reader with the given <code>Reader<code> for
	 * input.
	 *
	 * @param r		the <code>Reader</code> from which this network reader will
	 * 				parse a network
	 */
	public ExNewickReader(Reader r)
	{
		_reader = r;
	}

	/**
	 * This function reads a set of lines representing the network and builds
	 * the corresponding network.
	 *
	 * @return 						a network specified by the input from this
	 * 								network reader's <code>Reader</code>
	 * @throws IOException 			if the underlying <code>Reader</code>
	 * 								encounters an I/O error
	 * @throws ExNewickException 	if the input has the wrong format
	 */
	public Network<T> readNetwork() throws ExNewickException, IOException
	{
		// Read the input.
		getInput();

		// Parse the input to obtain the set of network nodes.
		_network_nodes = new LinkedList<BniNetNode<T>>();
		for (String str : _next_lines) {
			int i = str.indexOf('=');
			assert i >= 0 : "missing = character in subtree!";
			String name = str.substring(0, i);	// Get the name of the network node.
			BniNetNode<T> node = new BniNetNode<T>();

			name = name.trim();
			node.setName(name);
			_network_nodes.add(node);
		}

		// Build the network from the input.
		BniNetNode<T> root = parseTree(_first_line);	// Parse the main tree.
		for (String str : _next_lines) {
			parseTree(str);				// Parse subtrees.
		}

		// Check for node name uniqueness.
		_net = new BniNetwork<T>(root);
		if (_net.hasDuplicateNames()) {
			throw new ExNewickException("The network read contains duplicate names");
		}

		//Networks.removeBinaryNodes(_net);
		return _net;
	}


	/**
	 * Reads the input from the reader. Input strings are stored in
	 * _first_line and _next_lines.
	 *
	 * @throws ExNewickException	if the network stored in the reader is in
	 * 								some way malformed
	 *
	 * @throws IOException			if the underlying <code>Reader</code>
	 * 								encounters an I/O error
	 */
	private void getInput() throws IOException, ExNewickException
	{
		BufferedReader br = new BufferedReader(_reader);
        
		// Read the first non-blank line as the main tree.
		while ((_first_line = br.readLine()) != null) {
			_first_line = _first_line.trim();
			if (_first_line.length() != 0) {
				break;
			}
		}

		if (_first_line.indexOf(";") < 0) {
			_first_line += ";";
		}

		// Add name of the root.
		if (_first_line.indexOf("=") < 0) {
			String prefix = "N";
			int number = 0;
			String name = new String();

			do {
				name = prefix + (number++);
			} while (_first_line.contains(name));
			_first_line = name + "=" + _first_line;
		}

		// Read subsequent non-blank lines as subtrees.
		String str;
		_next_lines = new LinkedList<String>();

		while ((str = br.readLine()) != null) {
			str = str.trim();
			if (str.length() != 0) {
				if (str.indexOf("=") < 0) {
					throw new ExNewickException("= expected in subtree");
				}
				if (str.indexOf(";") < 0) {
					str += ";";
				}
				_next_lines.add(str);
			}
		}
	}

	/**
	 * Builds a tree from an input string.
	 *
	 * @param line	the string this method will parse.
	 *
	 * @return a BniNetNode which is the root of the tree.
	 *
	 * @throws ExNewickException	if the line contains unexpected characters
	 * 								or is otherwise malformed
	 * @throws IOException 			if the underlying <code>Reader</code>
	 * 								encounters an I/O error
	 */
	private BniNetNode<T> parseTree(String line) throws ExNewickException, IOException
	{
		// Set up the StreamTokenizer.
		StreamTokenizer tokenizer = new StreamTokenizer(new StringReader(line));
		initTokenizer(tokenizer);

		// Build the tree from line.
		BniNetNode<T> root = null;
		String root_name = null;
		Stack<BniNetNode<T>> nodes = new Stack<BniNetNode<T>>();

		// Get the name of the node which is the root of the tree contained in this line.
		if (tokenizer.nextToken() == StreamTokenizer.TT_WORD) {
			root_name = tokenizer.sval;
			if (tokenizer.nextToken() != '=') {
				throw new ExNewickException("= expected");
			}
		}
		else {
			throw new ExNewickException("Tree name expected");
		}

		// Initialize the nodes stack with the node corresponding to root_name.
		if ((root = getNetworkNode(root_name)) == null) {
			root = new BniNetNode<T>();
			root.setName(root_name);
			nodes.push(root);

		}
		//nodes.push(root);

		else{
			nodes.push(root);
			BniNetNode<T> nnode = new BniNetNode<T>();
			root.adoptChild(nnode, BniNetNode.NO_DISTANCE);
			nodes.push(nnode);
		}



		// Read leaves and internal nodes.
		double distance = BniNetNode.NO_DISTANCE;
		String name = BniNetNode.NO_NAME;
		double gamma = 1.0;
		int tt;
		BniNetNode<T> child;

		while ((tt = tokenizer.nextToken()) != StreamTokenizer.TT_EOF) {
			switch (tt) {
			case '(':
				nodes.push(new BniNetNode<T>());	// Create a spare node.
				break;
			case StreamTokenizer.TT_WORD:
				name = parseNodeName(tokenizer);
				if ((child = getNetworkNode(name)) != null) {	// This name is of a network node
					nodes.pop();
					nodes.push(child);	// Replace the spare node with that network node.
				}
				break;
			case ':':
				if ((tt = tokenizer.nextToken()) != StreamTokenizer.TT_WORD) {
					throw new ExNewickException("Number expected");
				}
				else {
					distance = parseNodeDistance(tokenizer);
				}
				break;
			case '^':
				if ((tt = tokenizer.nextToken()) != StreamTokenizer.TT_WORD) {
					throw new ExNewickException("Number expected");
				}
				else {
					gamma = parseNodeGamma(tokenizer);
				}
				break;
			case ',':
			case ')':
				// Join this node (as a child) to its parent.
				child = nodes.pop();
				child.setName(name);
				BniNetNode<T> parent = nodes.peek();
                parent.adoptChild(child, distance);
                child.setParentProbability(parent, gamma);

				// Prepare a spare node for the next child.
				if (tt == ',') {
					nodes.push(new BniNetNode<T>());
				}

				// Invalidate the old name and distance values.
				name = BniNetNode.NO_NAME;
				distance = BniNetNode.NO_DISTANCE;
				gamma = 1.0;
				break;
			case ';':
				if ((tt = tokenizer.nextToken()) != StreamTokenizer.TT_EOF) {
					throw new ExNewickException("Unexpected character after semicolon");
				}
				else {
					tokenizer.pushBack();
				}

				// The only element in the stack must be the tree's root.

				if(getNetworkNode(root_name)!= null){
					BniNetNode<T> nnode = nodes.pop();
					distance = Math.max(0, distance);
					nnode.setParentDistance(root, distance);
				}


				if (!nodes.pop().equals(root) || !nodes.isEmpty()) {
					throw new ExNewickException("Wrong eNewick format");
				}

				if (!name.equals(BniNetNode.NO_NAME)) {
					root.setName(name);
				}

				break;
			default:
				throw new ExNewickException("Unexpected character");
			}
		}

		return root;
	}

	/**
	 * Initializes a <code>StreamTokenizer</code> for parsing input.
	 *
	 * @param st 	a <code>StreamTokenizer</code>
	 */
	private void initTokenizer(StreamTokenizer st)
	{
		st.resetSyntax();

		// For parsing names.
		st.wordChars('A', 'Z');
		st.wordChars('a', 'z');
		st.wordChars('_', '_');

		// For parsing numbers.
		st.wordChars('0', '9');
		st.wordChars('-', '-');
		st.wordChars('.', '.');

		// Whitespaces.
		st.whitespaceChars(' ', ' ');
		st.whitespaceChars('\t', '\t');

		st.ordinaryChar('=');
		st.ordinaryChar('(');
		st.ordinaryChar(')');
		st.ordinaryChar(',');
		st.ordinaryChar(':');
		st.ordinaryChar(';');
		st.ordinaryChar('^');
	}

	/**
	 * Gets a network node with <code>name</code> stored in _network_nodes.
	 *
	 * @param name	the name of the network node.
	 * @return 		the corresponding network node in _network_nodes, or
	 * 				<code>null</code> if it is not found
	 */
	private BniNetNode<T> getNetworkNode(String name)
	{
		if (_network_nodes == null) {
			return null;	// Not found.
		}
		else {
			for (BniNetNode<T> node : _network_nodes) {
				if (node.getName().equals(name)) {
					return node;	// The node with the desired name is found.
				}
			}

			return null;	// Not found.
		}
	}

	/**
	 * Returns a node name from the supplied input.
	 *
	 * @param st	a <code>StreamTokenizer</code> to supply the input
	 * @return 		a name of the node
	 *
	 * @throws ExNewickException	if the last encountered token was not a
	 * 								valid node name character
	 * @throws IOException			if the supplied <code>StreamTokenizer</code>
	 * 								encounters an I/O error
	 */
	private String parseNodeName(StreamTokenizer st) throws ExNewickException, IOException
	{
		// Retrieve the name.
		String name = st.sval;

		// Enforce the extended newick format.
		int tt = st.nextToken();

		if (tt != ':' && tt != ',' && tt != ')' && tt != ';' && tt != '^') {
			throw new ExNewickException("Unexpected character");
		}
		else {
			st.pushBack();
		}

		return name;
	}

	/**
	 * Reads the distance for a node in the network.
	 *
	 * @param st 	a <code>StreamTokenizer</code> to supply the input
	 * @return		the distance value for a node.
	 *
	 * @throws ExNewickException	if a non-numerical or otherwise invalid
	 * 								character is encountered
	 * @throws IOException 			if the supplied <code>StreamTokenizer</code>
	 * 								encounters an I/O error
	 */
	private double parseNodeDistance(StreamTokenizer st) throws ExNewickException, IOException
	{
		// Retrieve the distance.
		double distance = BniNetNode.NO_DISTANCE;
		try {
			distance = Double.parseDouble(st.sval);
		}
		catch (NumberFormatException e) {
			throw new ExNewickException("Number expected"); // TODO Give the user more useful feedback about malformed networks (e.g. line numbers, etc.).
		}

		// Enforce the extended newick format.
		int tt = st.nextToken();

		if (tt != ',' && tt != ')' && tt != ';' && tt != '^') {
			throw new ExNewickException("Unexpected character");
		}
		else {
			st.pushBack();
		}

		return distance;
	}

	/**
	 * Reads the proportion for a network node in the network.
	 *
	 * @param st 	a <code>StreamTokenizer</code> to supply the input
	 * @return		the distance value for a node.
	 *
	 * @throws ExNewickException	if a non-numerical or otherwise invalid
	 * 								character is encountered
	 * @throws IOException 			if the supplied <code>StreamTokenizer</code>
	 * 								encounters an I/O error
	 */
	private double parseNodeGamma(StreamTokenizer st) throws ExNewickException, IOException
	{
		// Retrieve the distance.
		double gamma = 1.0;
		try {
			gamma = Double.parseDouble(st.sval);
		}
		catch (NumberFormatException e) {
			throw new ExNewickException("Number expected"); // TODO Give the user more useful feedback about malformed networks (e.g. line numbers, etc.).
		}

		// Enforce the extended newick format.
		int tt = st.nextToken();

		if (tt != ',' && tt != ')' && tt != ';') {
			throw new ExNewickException("Unexpected character");
		}
		else {
			st.pushBack();
		}

		return gamma;
	}

	// Data members
	private Reader _reader;		// Input stream.
	private Network<T> _net;	// The network read.
	private List<BniNetNode<T>> _network_nodes;	// Stores network nodes.
	private String _first_line;
	private List<String> _next_lines;
}

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

package edu.rice.cs.bioinfo.programs.phylonet.structs.sequence.model;

import java.io.*;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/31/11
 * Time: 3:32 PM
 * To change this template use File | Settings | File Templates.
 */
public class SequenceAlignment {

    /**
	 * This constructor constructs a sequence alignment from an array of taxon names
	 * and an array of strings representing sequences.
	 *
	 * NOTE: This function *assumes* that all strings in the array <code>sequences</code>
	 * have the same length.
	 *
	 * @param taxa 			an array of taxon names.
	 * @param sequences 	an array of strings representing sequences.
	 */
	public SequenceAlignment(String taxa[], String sequences[])
	{
		assert(taxa.length == sequences.length && taxa.length > 0);

		_taxa = taxa;
		_sequences = sequences;

		_size = _taxa.length;				// The number of taxa.
		_length = _sequences[0].length();	// Sequence length.
	}

	/**
	 * Instantiates a new, empty sequence alignment.
	 */
	public SequenceAlignment()
	{
		_size = _length = -1;
		_taxa = _sequences = null;
	}

	/**
	 * Parse the raw input sequence input. Input may be formatted in FASTA or
	 * PHYLIP format.
	 *
	 * @param  input				a set of taxa and sequences
	 *
	 * @throws SequenceException 	if the supplied sequences are improperly
	 * 								formatted
	 */
	public void readSequences(String input) throws SequenceException
	{
		// Check if we will read FASTA or PHYLIP format.
		boolean fasta = (input.charAt(0) == '>');
		StringReader sr = new StringReader(input);

		try {
			if (fasta) {
				readFastaSequences(sr);
			}
			else {
				readPhylipSequences(sr);
			}
		}
		catch (IOException e) {
			// we should never encounter an IOException, since the Reader is a simple StringReader
			throw new RuntimeException("supposedly impossible I/O error from a StringReader", e);
		}
		finally {
			sr.close();
		}
	}

	/**
	 * Implementation to read sequences in PHYLIP format (either sequential or interleaved).
	 * The PHYLIP format is as follows:
	 *
	 * <pre>
	 * #taxa	sequence-length	[I (interleaved)]
	 * taxa-name1	sequence1
	 * taxa-name2	sequence2
	 * ...
	 * 		sequence1'
	 * 		sequence2'
	 * 		... (if interleaved)</pre>
	 *
	 * In the PHYLIP format, the length of taxon names is 10, including blanks.
	 *
	 * @param r							the reader from which the sequences are
	 * 									to be read
	 *
	 * @throws IOException				if the supplied reader encounters an
	 * 									I/O error
	 * @throws SequenceException 		if the supplied sequences are
	 * 									improperly formatted
	 */
	public void readPhylipSequences(Reader r) throws IOException, SequenceException
	{
		BufferedReader br = new BufferedReader(r);
		boolean interleaved = false;

		try {

			// Read the header which contains the number of taxa and the length of DNA/protein sequences
			// and possibly I (for interleaved).
			String header = null;
			while ((header = br.readLine()) != null) {
				header = header.trim();
				if (header.length() != 0) {
					break;
				}
			}

			StreamTokenizer tokenizer = new StreamTokenizer(new StringReader(header));

			tokenizer.resetSyntax();
			tokenizer.wordChars('0', '9');
			tokenizer.whitespaceChars(' ', ' ');
			tokenizer.whitespaceChars('\t', '\t');

			try {
				int tt = tokenizer.nextToken();

				if (tt != StreamTokenizer.TT_WORD) {
					throw new SequenceException("Number of taxa expected");
				}
				else {
					try {
						_size = Integer.parseInt(tokenizer.sval);
						if (_size <= 0) {
							throw new SequenceException("A positive number expected");
						}
					}
					catch (NumberFormatException e) {
						throw new SequenceException("Number of taxa expected");
					}
				}

				tt = tokenizer.nextToken();
				if (tt != StreamTokenizer.TT_WORD) {
					throw new SequenceException("Sequence length expected");
				}
				else {
					try {
						_length = Integer.parseInt(tokenizer.sval);
						if (_length <= 0) {
							throw new SequenceException("A positive number expected");
						}
					}
					catch (NumberFormatException e) {
						throw new SequenceException("Sequence length expected");
					}
				}

				tt = tokenizer.nextToken();
				if (tt == StreamTokenizer.TT_WORD) {
					if (tokenizer.sval.equals("I")) {
						interleaved = true;
					}
					else {
						throw new SequenceException("Unexpected character in the header line of the sequence file");
					}

					tt = tokenizer.nextToken();
					if (tt != StreamTokenizer.TT_EOF) {
						throw new SequenceException("Unexpected character in the header line of the sequence file");
					}
				}
				else if (tt != StreamTokenizer.TT_EOF) {
					throw new SequenceException("Unexpected character in the header line of the sequence file");
				}
			}
			catch (IOException e) {
				// an I/O error here represents a bug, since we are just parsing our own string
				throw new RuntimeException("supposedly impossible I/O error!", e);
			}

			// Read sequences.
			StringBuffer lines[] = new StringBuffer[_size];
			int count = 0;

			String str = null;
			while ((str = br.readLine()) != null && count < _size) {
				str.trim();
				if (str.length() != 0) {	// Skip over all blank lines.
					lines[count] = new StringBuffer(str);
					count++;
				}
			}

			if (count < _size) {
				throw new SequenceException("The number of sequences does not match the number specified in the header line");
			}

			_taxa = new String[_size];

			for (int i = 0; i < _size; i++) {
				_taxa[i] = lines[i].substring(0, _phylip_name_length);
				_taxa[i] = _taxa[i].trim();

				lines[i].delete(0, _phylip_name_length);
			}

			// If interleaved = true, read next sequences.
			if (interleaved) {
				str = null;
				count = 0;
				while ((str = br.readLine()) != null) {
					str.trim();
					if (str.length() != 0) {
						lines[count].append(str);
						count++;
					}

					if (count == _size) {
						count = 0;
					}
				}
			}

			// Remove blank spaces from sequences.
			for (int i = 0; i < _size; i++) {
				for (int j = 0; j < lines[i].length(); j++) {
					if (lines[i].charAt(j) == ' ' || lines[i].charAt(j) == '\t' || lines[i].charAt(j) == '\n') {
						lines[i].deleteCharAt(j);
						j--;
					}
				}

				if (lines[i].length() != _length) {
					throw new SequenceException("The sequence length does not match the value specified in the header line");
				}
			}

			_sequences = new String[_size];
			for (int i = 0; i < _size; i++) {
				_sequences[i] = lines[i].toString();
			}

		}
		finally {
			// ensure reader closed in case we encountered an error
			br.close();
		}
	}

	/**
	 * Reads sequences, together with their corresponding taxa names.
	 *
	 * NOTE: The current format for the input is as follows:
	 * <pre>
	 * 	#taxa		sequence-length
	 * 	taxa-name1	sequence1
	 * 	taxa-name2	sequence2
	 * 	...</pre>
	 *
	 * @param r 	the reader from which this function gets the sequences
	 *
	 * @throws SequenceException	if the sequences or taxa are improperly
	 * 								formatted or otherwise malformed
	 * @throws IOException 			if the supplied <code>Reader</code>
	 * 								encounters an I/O error
	 */
	public void readPhyloNetSequences(Reader r) throws SequenceException, IOException
	{
		BufferedReader br = new BufferedReader(r);

		try {

			// Read the header which contains the number of taxa and sequence length.
			String header;
			while ((header = br.readLine()) != null) {
				header = header.trim();
				if (header.length() != 0) {
					break;
				}
			}

			// Get the number of taxa and the sequence length.
			StreamTokenizer tokenizer = new StreamTokenizer(new StringReader(header));

			tokenizer.resetSyntax();
			tokenizer.wordChars('0', '9');
			tokenizer.whitespaceChars(' ', ' ');
			tokenizer.whitespaceChars('\t', '\t');

			try {
				int tt = tokenizer.nextToken();

				if (tt != StreamTokenizer.TT_WORD) {
					throw new SequenceException("Number of taxa expected");
				}
				else {
					try {
						_size = Integer.parseInt(tokenizer.sval);
						if (_size <= 0) {
							throw new SequenceException("A positive number expected");
						}
					}
					catch (NumberFormatException e) {
						throw new SequenceException("Size expected");
					}
				}

				tt = tokenizer.nextToken();
				if (tt != StreamTokenizer.TT_WORD) {
					throw new SequenceException("Sequence length expected");
				}
				else {
					try {
						_length = Integer.parseInt(tokenizer.sval);
						if (_length <= 0) {
							throw new NumberFormatException("A positive number expected");
						}
					}
					catch (NumberFormatException e) {
						throw new SequenceException("Length expected");
					}
				}

				tt = tokenizer.nextToken();
				if (tt != StreamTokenizer.TT_EOF) {
					throw new SequenceException("Unexpected character");
				}
			}
			catch (IOException e) {
				// it's a bug if we encounter an I/O error here, since we're parsing our own string
				throw new RuntimeException("supposedly impossible I/O error parsing header string");
			}

			// Read the sequences.
			String lines[] = new String[_size];

			String str;
			int j = 0;

			while ((str = br.readLine()) != null) {
				str = str.trim();
				if (str.length() != 0) {
					if (j < _size) {
						lines[j] = str;
						j++;
					}
					else {
						throw new SequenceException("There are more taxa than the number specified in the header");
					}
				}
			}

			// Parse the sequences.
			_taxa = new String[_size];
			_sequences = new String[_size];

			for (int i = 0; i < _size; i++) {
				tokenizer = new StreamTokenizer(new StringReader(lines[i]));
				int tt = tokenizer.nextToken();

				if (tt != StreamTokenizer.TT_WORD) {
					throw new SequenceException("Taxon name expected");
				}
				else {
					_taxa[i] = tokenizer.sval;
				}

				tt = tokenizer.nextToken();
				if (tt != StreamTokenizer.TT_WORD) {
					throw new SequenceException("Sequence expected");
				}
				else if (tokenizer.sval.length() != _length) {
					throw new SequenceException("The length of the sequence does not match the length specified in the header");
				}
				else {
					_sequences[i] = tokenizer.sval;
				}

				tt = tokenizer.nextToken();
				if (tt != StreamTokenizer.TT_EOF) {
					throw new SequenceException("Unexpected character");
				}
			}

		}
		finally {
			// ensure closing in case we encountered an I/O error somewhere
			br.close();
		}
	}

	/**
	 * Reads a set of sequences in the FASTA format from the supplied
	 * <code>Reader</code>.
	 *
	 * @param 	r	the <code>Reader</code> from which the sequences will be
	 * 				read
	 *
	 * @throws SequenceException if the sequences are in some manner malformed
	 * @throws IOException 		 if the supplied <code>Reader</code> encounters
	 * 							 an I/O error
	 */
	public void readFastaSequences(Reader r) throws SequenceException, IOException
	{
		// There are no header in the FASTA format. Instead, each line begins with a '>' marker, followed
		// by a sequence identifier.
		BufferedReader br = new BufferedReader(r);
		List<StringBuffer> lines = new LinkedList<StringBuffer>();
		List<String> seqIds = new LinkedList<String>();
		String str = null;
		StringBuffer seq = null;

		try {
			while ((str = br.readLine()) != null) {
				str.trim();
				if (str.length() != 0) {
					if (str.charAt(0) == '>') {	// Begin of a new sequence.
						int i = str.indexOf('|');
						if (i < 0) {
							i = str.length();
						}
						seqIds.add(str.substring(1, i));

						if (seq != null) {
							lines.add(seq);	// Record
						}
						seq = new StringBuffer();
					}
					else {	// Assemble a new sequence.
						seq.append(str);
					}
				}
			}
			lines.add(seq);	// Append the last sequence.
		}
		finally {
			br.close();
		}

		// Check the sequence length.
		_size = lines.size();
		if (_size <= 0) {
			throw new SequenceException("There are no sequences in the file");
		}

		_taxa = new String[_size];
		_sequences = new String[_size];

		_length = lines.get(0).length();
		for (int i = 0; i < _size; i++) {
			if (lines.get(i).length() != _length) {
				throw new SequenceException("Sequences have different lengths");
			}
			else {
				_taxa[i] = seqIds.get(i);
				_sequences[i] = lines.get(i).toString();
			}
		}

	}

	/**
	 * Returns the number of sequences.
	 *
	 *  @return the number of sequences
	 */
	public int getSequenceCount()
	{
		return _size;
	}

	/**
	 * Returns the sequence length.
	 *
	 * @return the sequence length
	 */
	public int getSequenceLength()
	{
		return _length;
	}

	/**
	 * Returns the taxon names of this sequence alignment.
	 *
	 * @return the taxon names of this sequence alignment
	 */
	public String[] getTaxa()
	{
		return _taxa;
	}

	/**
	 * Returns the sequences contained by this sequence alignment.
	 *
	 * @return	the sequences contained by this sequence alignment
	 */
	public String[] getSequences()
	{
		return _sequences;
	}

	/**
	 * Returns an array of blocks extracted from _sequences. The block extends from start to stop - 1.
	 * Therefore, the block size is stop - start.
	 *
	 * NOTE: To get the entire sequence, assign values to start and stop such that stop < start.
	 *
	 * @param start 	indicates the start of the block.
	 * @param stop 		indicates the stop of the block.
	 *
	 * @return	the specified subset of sequences
	 */
	public String[] getBlock(int start, int stop)
	{
		assert(_size > 0 && _length > 0);
		assert(_taxa != null && _sequences != null);

		if (stop < start) {
			return _sequences;
		}
		else if (start < 0 || start >= _length || stop < 0 || stop > _length) {
			throw new IllegalArgumentException("Wrong index values for the start and stop of a block.");
		}
		else {
			String temp[] = new String[_size];

			for (int i = 0; i < _size; i++) {
				temp[i] = _sequences[i].substring(start, stop);
			}
			return temp;
		}
	}

	// Data members
	private int _length;		// The length of the sequences in this sequence alignment.
	private int _size;			// The number of sequences in this alignment. It is also the number of taxa.

	private String _taxa[];	// Holds the taxon names.
	private String _sequences[];	// Holds sequences as string of characters.
	private static int _phylip_name_length = 10;

}

package edu.rice.bioinfo.programs.richnewickwebview;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.table.AbstractTableModel;

import edu.rice.bioinfo.library.language.richnewick._1_0.RichNewickReadError;
import edu.rice.bioinfo.library.language.richnewick._1_0.RichNewickReadException;
import edu.rice.bioinfo.library.language.richnewick._1_0.ast.DAGFactory;
import edu.rice.bioinfo.library.language.richnewick._1_0.ast.Network;
import edu.rice.bioinfo.library.language.richnewick._1_0.csa.CSAError;
import edu.rice.bioinfo.library.language.richnewick._1_0.graphbuilding.GraphBuilderDOT;
import edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast.*;

import com.mxgraph.swing.mxGraphComponent;
import com.mxgraph.view.mxGraph;
import edu.rice.bioinfo.library.programming.Func;
import edu.rice.bioinfo.library.programming.Func1;
import edu.uci.ics.jung.algorithms.layout.*;
import edu.uci.ics.jung.algorithms.layout.SpringLayout;
import edu.uci.ics.jung.graph.DirectedGraph;
import sun.org.mozilla.javascript.internal.Function;

import java.applet.Applet;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;

public class Program extends JApplet
{

    /**
     *
     */
    private static final long serialVersionUID = -2707712944901661771L;

    private Network _currentNetwork;

    @Override
    public void start()
    {
        try
        {



            JPanel inputPane = new JPanel();
            inputPane.setLayout(new BorderLayout());



            final JTextArea textInput = new JTextArea();
            textInput.setEditable(true);
            textInput.setText("((1, ((2, (3, (4)Y#H1)g)e, (((Y#H1, 5)h, 6)f)X#H2)c)a, ((X#H2, 7)d, 8)b)r;");
            inputPane.add(textInput, BorderLayout.CENTER);

            JPanel nonTextInput = new JPanel();
            nonTextInput.setLayout(new BorderLayout());
            inputPane.add(nonTextInput, BorderLayout.SOUTH);

            final JComboBox layouts = new JComboBox(new Object[]
                    {
                            new Func1<DirectedGraph<Object,Object>, Layout<Object,Object>>()
                            {
                                public Layout<Object,Object> execute(DirectedGraph<Object,Object> dag) {
                                    return new CircleLayout<Object, Object>(dag);
                                }

                                @Override public String toString()
                                {
                                    return "Circle Layout";
                                }
                            },


                            new Func1<DirectedGraph<Object,Object>, Layout<Object,Object>>()
                            {
                                public Layout<Object,Object> execute(DirectedGraph<Object,Object> dag) {
                                    return new DAGLayout<Object,Object>(dag);
                                }

                                @Override public String toString()
                                {
                                    return "DAG Layout";
                                }
                            },

                            new Func1<DirectedGraph<Object,Object>, Layout<Object,Object>>()
                            {
                                public Layout<Object,Object> execute(DirectedGraph<Object,Object> dag) {
                                    return new FRLayout<Object, Object>(dag);
                                }

                                @Override public String toString()
                                {
                                    return "FR Layout";
                                }
                            },

                            new Func1<DirectedGraph<Object,Object>, Layout<Object,Object>>()
                            {
                                public Layout<Object,Object> execute(DirectedGraph<Object,Object> dag) {
                                    return new FRLayout2<Object, Object>(dag);
                                }

                                @Override public String toString()
                                {
                                    return "FR 2 Layout";
                                }
                            },

                            new Func1<DirectedGraph<Object,Object>, Layout<Object,Object>>()
                            {
                                public Layout<Object,Object> execute(DirectedGraph<Object,Object> dag) {
                                    return new ISOMLayout<Object, Object>(dag);
                                }

                                @Override public String toString()
                                {
                                    return "ISOM Layout";
                                }
                            },

                            new Func1<DirectedGraph<Object,Object>, Layout<Object,Object>>()
                            {
                                public Layout<Object,Object> execute(DirectedGraph<Object,Object> dag) {
                                    return new KKLayout<Object, Object>(dag);
                                }

                                @Override public String toString()
                                {
                                    return "KK Layout";
                                }
                            },

                            new Func1<DirectedGraph<Object,Object>, Layout<Object,Object>>()
                            {
                                public Layout<Object,Object> execute(DirectedGraph<Object,Object> dag) {
                                    return new SpringLayout<Object, Object>(dag);
                                }

                                @Override public String toString()
                                {
                                    return "Spring Layout";
                                }
                            },

                            new Func1<DirectedGraph<Object,Object>, Layout<Object,Object>>()
                            {
                                public Layout<Object,Object> execute(DirectedGraph<Object,Object> dag) {
                                    return new SpringLayout2<Object, Object>(dag);
                                }

                                @Override public String toString()
                                {
                                    return "Spring 2 Layout";
                                }
                            }


                    });
            nonTextInput.add(layouts, BorderLayout.CENTER);

            JPanel buttonPanel = new JPanel();
            buttonPanel.setLayout(new BorderLayout());
            nonTextInput.add(buttonPanel, BorderLayout.EAST);

            JButton parse = new JButton("Parse");
            buttonPanel.add(parse, BorderLayout.WEST);

            final JButton export = new JButton("Show DOT");
            export.setEnabled(false);
            textInput.getDocument().addDocumentListener(new DocumentListener() {
                public void insertUpdate(DocumentEvent e) {
                    export.setEnabled(false);
                }

                public void removeUpdate(DocumentEvent e) {
                    export.setEnabled(false);
                }

                public void changedUpdate(DocumentEvent e) {
                    export.setEnabled(false);
                }
            });
            buttonPanel.add(export, BorderLayout.EAST);

            final JSplitPane mainPanel = new JSplitPane(JSplitPane.VERTICAL_SPLIT, inputPane, new JPanel());
            setContentPane(mainPanel);


            parse.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {

                    try
                    {
                        // mainPanel.setBottomComponent(new JPanel());
                        String networkText = textInput.getText();
                        JungGraphBuilder builder = new JungGraphBuilder( (Func1<DirectedGraph<Object,Object>,Layout<Object,Object>>) layouts.getSelectedItem());


                        RichNewickReaderAST_ANTLR reader = new RichNewickReaderAST_ANTLR();

                        _currentNetwork = reader.read(new ByteArrayInputStream(networkText.getBytes()), builder);
                        mainPanel.setBottomComponent(builder.getVisualization());
                        export.setEnabled(true);
                    }
                    catch(final RichNewickReadException ex)
                    {
                        JPanel errorPanel = new JPanel();

                        JTable errorText = new JTable(new AbstractTableModel() {
                            public String getColumnName(int col) {

                                switch (col)
                                {
                                    case 0:
                                        return "Message";
                                    case 1:
                                        return "Line";
                                    case 2:
                                        return "Column";
                                    default:
                                        throw new IllegalArgumentException();

                                }
                            }
                            public int getRowCount() { return ex.Errors.size(); }
                            public int getColumnCount() { return 3; }
                            public Object getValueAt(int row, int col) {

                                RichNewickReadError error = ex.Errors.get(row);

                                 switch (col)
                                {
                                    case 0:
                                        return error.Message;
                                    case 1:
                                        return error.LineNumber;
                                    case 2:
                                        return error.ColumnNumber;
                                    default:
                                        throw new IllegalArgumentException();

                                }

                            }

                            public boolean isCellEditable(int row, int col)
                            { return false; }

                            public void setValueAt(Object value, int row, int col) {

                            }
                        });
                        errorText.setCellSelectionEnabled(false);

                        errorPanel.setLayout(new BorderLayout());
                        errorPanel.add(errorText.getTableHeader(), BorderLayout.NORTH);
                        errorPanel.add(errorText, BorderLayout.CENTER);

                        mainPanel.setBottomComponent(errorPanel);
                        export.setEnabled(false);
                        return;
                    }
                    catch(Exception ex)
                    {
                        JTextArea errorText = new JTextArea();
                        String errorMessage = ex.getMessage() == null ? ex.getClass().getName() : ex.getMessage();

                        for(StackTraceElement el: ex.getStackTrace())
                        {
                            errorMessage+= el.toString() + "\n";
                        }
                        errorText.setText(errorMessage);
                        errorText.setEditable(false);
                        mainPanel.setBottomComponent(errorText);
                        export.setEnabled(false);
                        return;
                    }



                }
            });

            export.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {

                    GraphBuilderDOT dotBuilder = new GraphBuilderDOT();
                    DAGFactory.makeDAG(_currentNetwork, dotBuilder);

                    JTextArea dotText = new JTextArea();
                    dotText.setEditable(false);
                    dotText.setText(dotBuilder.getDOT());
                    mainPanel.setBottomComponent(dotText);


                }
            });








        }
        catch(Exception e)
        {
            throw new RuntimeException(e);
        }



    }



}

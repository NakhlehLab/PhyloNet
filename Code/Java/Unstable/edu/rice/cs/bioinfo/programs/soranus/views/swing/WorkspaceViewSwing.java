package edu.rice.cs.bioinfo.programs.soranus.views.swing;

import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.*;
import edu.rice.cs.bioinfo.programs.soranus.views.WorkspaceView;

import javax.swing.*;
import javax.swing.tree.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.File;
import java.util.*;


/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/30/13
 * Time: 6:09 PM
 * To change this template use File | Settings | File Templates.
 */
public class WorkspaceViewSwing<DR extends DataRecord,
        AR extends AnalysisRecord> extends JFrame implements WorkspaceView<DR,AR>
{
    class RecordTreeNode<R> extends DefaultMutableTreeNode
    {
        public final R Record;

        public RecordTreeNode(String title, R record)
        {
            super(title);
            Record = record;
        }
    }

    private WorkspaceVM _viewModel;

    private JSplitPane _splitPane;

    private JTree _projectTree;

    private JPanel _documentPanel = new JPanel();

    private DefaultTreeModel _projectTreeModel;

    private DefaultMutableTreeNode _projectFolder = new DefaultMutableTreeNode("Project");

    private DefaultMutableTreeNode _dataFolder = new DefaultMutableTreeNode("Data");

    private DefaultMutableTreeNode _analysisFolder = new DefaultMutableTreeNode("Analysis");

    private Set<TreeNode> _dataRecordTreeNodes = new HashSet<TreeNode>();

    public final Observable1<File> AddDataRequested = new Observable1<File>();

    public Observable1<File> getAddDataRequested()
    {
        return AddDataRequested;
    }

    public final Observable1<DR> NeighborJoiningAnalysisRequested = new Observable1<DR>();

    public Observable1<DR> getNeighborJoiningAnalysisRequested()
    {
        return NeighborJoiningAnalysisRequested;
    }

    public final Observable3<DR,DR,DR> SnitkinTransMapAnalysisRequested = new Observable3<DR,DR,DR>();

    public Observable3<DR,DR,DR> getSnitkinTransMapAnalysisRequested() {
        return SnitkinTransMapAnalysisRequested;
    }

    public final Observable1<DR> DetectRecombAnalysisRequested = new Observable1<DR>();

    public Observable1<DR> getDetectRecombAnalysisRequested()
    {
        return DetectRecombAnalysisRequested;
    }

    private Map<DefaultMutableTreeNode,AR> _analysisTreeNodeToRecord = new HashMap<DefaultMutableTreeNode, AR>();

    private Map<DefaultMutableTreeNode,DR> _dataTreeNodeToRecord = new HashMap<DefaultMutableTreeNode, DR>();

    public final Observable1<DR> MinSpanTreeSnpRequested = new Observable1<DR>();

    public Observable1<DR> getMinSpanTreeSnpRequested()
    {
        return MinSpanTreeSnpRequested;
    }

    public final Observable1<AR> AnalysisRecordSelected = new Observable1<AR>();

    public Observable1<AR> getAnalysisRecordSelected()
    {
        return AnalysisRecordSelected;
    }

    public final Observable1<DR> DataRecordSelected = new Observable1<DR>();

    public Observable1<DR> getDataRecordSelected()
    {
        return DataRecordSelected;
    }

    public void startView() {
        this.setExtendedState(this.getExtendedState() | JFrame.MAXIMIZED_BOTH);
        this.setVisible(true);
    }


    private JButton _createSequencingsData = new JButton("Create Sequencings from VAAL");

    public WorkspaceViewSwing(WorkspaceVM viewModel)
    {
        _viewModel= viewModel;

        _viewModel.addDataRecordAddedListener(new Proc1<DR>() {
            public void execute(DR record) {
                onDataRecordAdded(record);
            }
        });
        _viewModel.addFocusDocumentChangedListener(new Proc1<DocumentVM>() {
            public void execute(DocumentVM documentVM) {
                onFocusDocumentChanged(documentVM);
            }
        });
        _viewModel.addAnalysisAddedListener(new Proc1<AR>() {
            public void execute(AR record) {
                onAnalysisAdded(record);
            }
        });

        this.setTitle(viewModel.Title);
        this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        _projectFolder.add(_dataFolder);
        _projectFolder.add(_analysisFolder);


        _projectTreeModel = new DefaultTreeModel(_projectFolder);
        _projectTree = makeProjectTree(_projectTreeModel);
        //   _projectTreeModel.reload();

        _documentPanel.setLayout(new BorderLayout());

        _splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, new JScrollPane(_projectTree), _documentPanel);
        _splitPane.resetToPreferredSizes();
        this.getContentPane().add(_splitPane);

        JToolBar menuBar = new JToolBar();

        menuBar.add(_createSequencingsData);
        this.getContentPane().add(menuBar, BorderLayout.NORTH);


    }

    public void addCreateSequencingsDataRequestListener(final Proc listener)
    {
        _createSequencingsData.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                listener.execute();;
            }
        });
    }



    private void onFocusDocumentChanged(DocumentVM documentVM)
    {
        JComponent newFocusView = documentVM.execute(new DocumentVMAlgo<JComponent, RuntimeException>()
        {
            public <N, E> JPanel forTransMapVM(TransMapVM<N,E> vm) throws RuntimeException {
                return new TransMapViewDot<N,E>(vm);
            }

            public <N, E> JPanel forNeighborJoiningVM(NeighborJoiningVM<N,E> vm) throws RuntimeException {
                return new NeighborJoiningViewDot<N,E>(vm);
            }

            public <N, Ed> JComponent forTreeVM(TreeVM<N, Ed> vm) throws RuntimeException
            {
                return new TreeViewDot<N,Ed>(vm);
            }

            public JComponent forXMLDataVM(XMLDataVM xmlDataVM) {
                return new XMLDataView(xmlDataVM);
            }

            public JComponent forVAALOutDataVM(VAALOutDataVM vaalOutDataVM) {
                return new VAALOutDataView(vaalOutDataVM);
            }

            public JComponent forRecombResult(RecombResultVM vm) throws RuntimeException
            {
                return new RecombResultView(vm);
            }


        });
        _documentPanel.removeAll();
        _documentPanel.add(newFocusView, BorderLayout.CENTER);
        _documentPanel.validate();

    }

    private void onAnalysisAdded(AR record)
    {
        DefaultMutableTreeNode newEntryNode = new DefaultMutableTreeNode(record.Title);
        _analysisTreeNodeToRecord.put(newEntryNode, record);

        _analysisFolder.add(newEntryNode);
        _projectTreeModel.reload(_analysisFolder);

        TreePath pathFromRootToAnalysisFolder = getPathFromRootToFolder(_analysisFolder);
        _projectTree.expandPath(pathFromRootToAnalysisFolder);
        _splitPane.resetToPreferredSizes();
    }

    private TreePath getPathFromRootToFolder(DefaultMutableTreeNode folder)
    {
        LinkedList<TreeNode> pathFolderToRootAccum = new LinkedList<TreeNode>();
        pathFolderToRootAccum.add(folder);

        while(pathFolderToRootAccum.getLast() != _projectFolder)
        {
            pathFolderToRootAccum.add(pathFolderToRootAccum.getLast().getParent());
        }

        TreePath pathFromRootToFolder = new TreePath(pathFolderToRootAccum.removeLast()) ;

        while(!pathFolderToRootAccum.isEmpty())
        {
            pathFromRootToFolder = pathFromRootToFolder.pathByAddingChild(pathFolderToRootAccum.removeLast());
        }


        return pathFromRootToFolder;
    }


    private JTree makeProjectTree(TreeModel treeModel)
    {


        final JTree projectTree = new JTree(treeModel);





        DefaultTreeCellRenderer cellRenderer = new DefaultTreeCellRenderer()
        {
            @Override
            public Component getTreeCellRendererComponent(JTree tree, Object value, boolean sel, boolean expanded,
                                                          boolean leaf, int row, boolean hasFocus) {

                super.getTreeCellRendererComponent(
                        tree, value, sel,
                        expanded, leaf, row,
                        hasFocus);

                if(value == _dataFolder || value == _analysisFolder)
                {
                    setIcon(this.getDefaultOpenIcon());
                }

                return this;
            }
        };


        projectTree.getSelectionModel().setSelectionMode(TreeSelectionModel.CONTIGUOUS_TREE_SELECTION);
        projectTree.setCellRenderer(cellRenderer);
        projectTree.addMouseListener(new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {

                if(e.getSource() != projectTree)
                    return;

                TreePath path = projectTree.getPathForLocation ( e.getX (), e.getY () );
                Rectangle pathBounds = projectTree.getUI ().getPathBounds ( projectTree, path );

                if(pathBounds == null)
                    return;

                Object lastPathComponent = path.getLastPathComponent();

                if(SwingUtilities.isRightMouseButton(e))
                {
                    if(lastPathComponent == _dataFolder)
                    {
                        onDataFolderRightClick(e);
                    }
                    else if(_dataRecordTreeNodes.contains(lastPathComponent))
                    {
                        onDataRecordRightClick(e);
                    }
                }
                else if(SwingUtilities.isLeftMouseButton(e))
                {
                    if(_analysisTreeNodeToRecord.containsKey(lastPathComponent))
                    {
                        AR record = _analysisTreeNodeToRecord.get(lastPathComponent);
                        getAnalysisRecordSelected().notify(record);
                    }
                    else if(_dataTreeNodeToRecord.containsKey(lastPathComponent))
                    {
                        DR record = _dataTreeNodeToRecord.get(lastPathComponent);
                        getDataRecordSelected().notify(record);
                    }
                }
            }
        });

        return projectTree;
    }




    private void onDataRecordRightClick(MouseEvent e)
    {
        TreePath[] selections = _projectTree.getSelectionModel().getSelectionPaths();

        JPopupMenu analysisOptions = null;

        if(selections.length == 1)
        {
            final DR firstEntry =  (DR) ((RecordTreeNode<DR>)selections[0].getLastPathComponent()).Record;
            final boolean firstEntryTitleIsSequencingsData = _viewModel.isDataRecordTitleRepresentingSequencingsData(firstEntry);

            if(firstEntryTitleIsSequencingsData)
            {
                analysisOptions = new JPopupMenu ();
                final JMenuItem doNeighborJoin = new JMenuItem ( "Infer Neighbor Joining Tree" );
                doNeighborJoin.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        if(e.getSource() == doNeighborJoin)
                        {

                            getNeighborJoiningAnalysisRequested().notify(firstEntry);
                        }
                    }
                });
                analysisOptions.add (doNeighborJoin);


                final JMenuItem detectRecomb = new JMenuItem ( "Detect recombination" );
                detectRecomb.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        if(e.getSource() == detectRecomb)
                        {
                            getDetectRecombAnalysisRequested().notify(firstEntry);
                        }
                    }
                });
                analysisOptions.add (detectRecomb);

                final JMenuItem doMinimumSpanningTreeSnp =
                        new JMenuItem("Infer Min Span Tree (SNP)");
                doMinimumSpanningTreeSnp.addActionListener(new ActionListener()
                {
                    public void actionPerformed(ActionEvent e)
                    {
                        if(e.getSource() == doMinimumSpanningTreeSnp)
                        {
                            getMinSpanTreeSnpRequested().notify(firstEntry);

                        }
                    }
                });
                analysisOptions.add (doMinimumSpanningTreeSnp);


            }

        }
        else if(selections.length == 3)
        {
            final DR firstEntry =  (DR) ((RecordTreeNode<DR>)selections[0].getLastPathComponent()).Record;
            final boolean firstEntryTitleIsSequencingsData = _viewModel.isDataRecordTitleRepresentingSequencingsData(firstEntry);
            boolean firstEntryTitleIsTraceData = _viewModel.isDataRecordTitleRepresentingTraceData(firstEntry);

            final DR secondEntry = (DR) ((RecordTreeNode<DR>)selections[1].getLastPathComponent()).Record;
            boolean secondEntryTitleIsSequencingsData = _viewModel.isDataRecordTitleRepresentingSequencingsData(secondEntry);
            boolean secondEntryTitleIsTraceData = _viewModel.isDataRecordTitleRepresentingTraceData(secondEntry);

            final DR thirdEntry = (DR) ((RecordTreeNode<DR>)selections[2].getLastPathComponent()).Record;

            DR sequencingsEntry = null, traceEntry = null, firstPositiveEntry = null;
            for(DR entry : Arrays.asList(firstEntry, secondEntry, thirdEntry))
            {
                if(_viewModel.isDataRecordTitleRepresentingSequencingsData(entry))
                {
                    sequencingsEntry = entry;
                }
                else if(_viewModel.isDataRecordTitleRepresentingTraceData(entry))
                {
                    traceEntry = entry;
                }
                else if(_viewModel.isDataRecordTitleRepresentingFirstPositiveData(entry))
                {
                    firstPositiveEntry = entry;
                }
            }


            boolean showSnitkinAnalysisOption = sequencingsEntry != null && traceEntry != null &&
                    firstPositiveEntry != null;




            if(showSnitkinAnalysisOption)
            {
                final DR sequencingsEntryFinal = sequencingsEntry, traceEntryFinal = traceEntry,
                        firstPositiveEntryFinal = firstPositiveEntry;
                analysisOptions = new JPopupMenu ();
                final JMenuItem doSnitkin = new JMenuItem ( "Infer Snitkin Trans Map" );
                doSnitkin.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        if(e.getSource() == doSnitkin)
                        {
                            getSnitkinTransMapAnalysisRequested().notify(sequencingsEntryFinal, traceEntryFinal, firstPositiveEntryFinal);
                        }
                    }
                });
                analysisOptions.add (doSnitkin);
            }

        }

        if(analysisOptions != null)
        {
            analysisOptions.show (_projectTree, e.getX(), e.getY() );
        }
    }

    private void onDataFolderRightClick(MouseEvent e)
    {
        JPopupMenu menu = new JPopupMenu ();
        final JMenuItem addData = new JMenuItem ( "Add data..." );
        addData.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if(e.getSource() == addData)
                {
                    addDataRequested();
                }
            }
        });
        menu.add (addData);
        menu.show (_projectTree, e.getX(), e.getY() );
    }

    private void onDataRecordAdded(DR record)
    {
        RecordTreeNode<DR> newEntryNode = new RecordTreeNode<DR>(record.Title, record);
        _dataTreeNodeToRecord.put(newEntryNode, record);
        _dataFolder.add(newEntryNode);
        _dataRecordTreeNodes.add(newEntryNode);

        LinkedList<TreeNode> pathFromNewEntryToRootAccum = new LinkedList<TreeNode>();
        pathFromNewEntryToRootAccum.add(newEntryNode);

        while(pathFromNewEntryToRootAccum.getLast() != _projectFolder)
        {
            pathFromNewEntryToRootAccum.add(pathFromNewEntryToRootAccum.getLast().getParent());
        }

        pathFromNewEntryToRootAccum.removeFirst();

        TreePath pathFromRootToNewEntryParent = new TreePath(pathFromNewEntryToRootAccum.removeLast());

        while(!pathFromNewEntryToRootAccum.isEmpty())
        {
            pathFromRootToNewEntryParent = pathFromRootToNewEntryParent.pathByAddingChild(pathFromNewEntryToRootAccum.removeLast());
        }

        _projectTreeModel.reload(_dataFolder);
        _projectTree.expandPath(pathFromRootToNewEntryParent);
        _splitPane.resetToPreferredSizes();



    }

    private void addDataRequested()
    {

        JFileChooser fc = new JFileChooser();
        fc.setMultiSelectionEnabled(true);
        int returnVal = fc.showOpenDialog(this);

        if (returnVal == JFileChooser.APPROVE_OPTION) {
            for(File file : fc.getSelectedFiles())
            {
                AddDataRequested.notify(file);
            }

        }
    }

}

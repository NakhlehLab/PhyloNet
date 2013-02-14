package edu.rice.cs.bioinfo.programs.soranus.views.swing;

import edu.rice.cs.bioinfo.library.programming.Proc1;
import edu.rice.cs.bioinfo.library.programming.Proc3;
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
public class WorkspaceViewSwing extends JFrame implements WorkspaceView
{
    private WorkspaceVM _viewModel;

    private JSplitPane _splitPane;

    private JTree _projectTree;

    private JPanel _documentPanel = new JPanel();

    private DefaultTreeModel _projectTreeModel;

    private DefaultMutableTreeNode _projectFolder = new DefaultMutableTreeNode("Project");

    private DefaultMutableTreeNode _dataFolder = new DefaultMutableTreeNode("Data");

    private DefaultMutableTreeNode _analysisFolder = new DefaultMutableTreeNode("Analysis");

    private Set<TreeNode> _dataRecordTreeNodes = new HashSet<TreeNode>();

    private Set<Proc1<File>> _addDataListeners = new HashSet<Proc1<File>>();

    private Map<DefaultMutableTreeNode,WorkspaceVM.AnalysisRecord> _analysisTreeNodeToRecord = new HashMap<DefaultMutableTreeNode, WorkspaceVM.AnalysisRecord>();

    private Map<DefaultMutableTreeNode,WorkspaceVM.DataRecord> _dataTreeNodeToRecord = new HashMap<DefaultMutableTreeNode, WorkspaceVM.DataRecord>();

    public void addDataAddRequestListener(Proc1<File> listener)
    {
        _addDataListeners.add(listener);
    }

    private Set<Proc3<String,String,String>> _snitkinTransMapAnalysisRequestedListeners = new HashSet<Proc3<String,String,String>>();

    public void addSnitkinTransMapAnalysisRequestedListener(Proc3<String,String,String> listener)
    {
        _snitkinTransMapAnalysisRequestedListeners.add(listener);
    }

    private Set<Proc1<String>> _neighborJoiningAnalysisRequestedListeners = new HashSet<Proc1<String>>();

    public void addNeighborJoiningAnalysisRequestedListener(Proc1<String> listener)
    {
        _neighborJoiningAnalysisRequestedListeners.add(listener);
    }

    private Set<Proc1<WorkspaceVM.AnalysisRecord>> _analysisRecordSelectedListeners = new HashSet<Proc1<WorkspaceVM.AnalysisRecord>>();

    public void addAnalysisRecordSelectedListener(Proc1<WorkspaceVM.AnalysisRecord> listener)
    {
        _analysisRecordSelectedListeners.add(listener);
    }

    private Set<Proc1<WorkspaceVM.DataRecord>> _dataRecordSelectedListeners = new HashSet<Proc1<WorkspaceVM.DataRecord>>();

    public void addDataRecordSelectedListener(Proc1<WorkspaceVM.DataRecord> listener)
    {
        _dataRecordSelectedListeners.add(listener);
    }

    public void startView() {
        this.setExtendedState(this.getExtendedState() | JFrame.MAXIMIZED_BOTH);
        this.setVisible(true);
    }


    public WorkspaceViewSwing(WorkspaceVM viewModel)
    {
        _viewModel= viewModel;

        _viewModel.addDataRecordAddedListener(new Proc1<WorkspaceVM.DataRecord>() {
            public void execute(WorkspaceVM.DataRecord record) {
                onDataRecordAdded(record);
            }
        });
        _viewModel.addFocusDocumentChangedListener(new Proc1<DocumentVM>() {
            public void execute(DocumentVM documentVM) {
                onFocusDocumentChanged(documentVM);
            }
        });
        _viewModel.addAnalysisAddedListener(new Proc1<WorkspaceVM.AnalysisRecord>() {
            public void execute(WorkspaceVM.AnalysisRecord record) {
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

    }



    private void onFocusDocumentChanged(DocumentVM documentVM)
    {
        JComponent newFocusView = documentVM.execute(new DocumentVMAlgo<JComponent, RuntimeException>()
        {
            public <N, E> JPanel forTransMapVM(TransMapVM<N,E> vm) throws RuntimeException {
                return new TransMapViewJFiles<N,E>(vm);
            }

            public <N, E> JPanel forNeighborJoiningVM(NeighborJoiningVM<N,E> vm) throws RuntimeException {
                return new NeighborJoiningViewJFiles<N,E>(vm);
            }

            public JComponent forXMLDataVM(XMLDataVM xmlDataVM) {
                return new XMLDataView(xmlDataVM);
            }
        });
        _documentPanel.removeAll();
        _documentPanel.add(newFocusView, BorderLayout.CENTER);
        _documentPanel.validate();

    }

    private void onAnalysisAdded(WorkspaceVM.AnalysisRecord record)
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
                        WorkspaceVM.AnalysisRecord record = _analysisTreeNodeToRecord.get(lastPathComponent);
                        onAnalysisRecordLeftClick(record);
                    }
                    else if(_dataTreeNodeToRecord.containsKey(lastPathComponent))
                    {
                        WorkspaceVM.DataRecord record = _dataTreeNodeToRecord.get(lastPathComponent);
                        onDataRecordLeftClick(record);
                    }
                }
            }
        });

        return projectTree;
    }

    private void onAnalysisRecordLeftClick(WorkspaceVM.AnalysisRecord record) {

        for(Proc1<WorkspaceVM.AnalysisRecord> listener : _analysisRecordSelectedListeners)
        {
            listener.execute(record);
        }
    }

    private void onDataRecordLeftClick(WorkspaceVM.DataRecord record) {

        for(Proc1<WorkspaceVM.DataRecord> listener : _dataRecordSelectedListeners)
        {
            listener.execute(record);
        }
    }

    private void onDataRecordRightClick(MouseEvent e)
    {
        TreePath[] selections = _projectTree.getSelectionModel().getSelectionPaths();

        JPopupMenu analysisOptions = null;

        if(selections.length == 1)
        {
            final String firstEntryTitle =  (String) ((DefaultMutableTreeNode)selections[0].getLastPathComponent()).getUserObject();
            final boolean firstEntryTitleIsSequencingsData = _viewModel.isDataRecordTitleRepresentingSequencingsData(firstEntryTitle);

            if(firstEntryTitleIsSequencingsData)
            {
                final String sequencingsTitleFinal = firstEntryTitle;
                analysisOptions = new JPopupMenu ();
                final JMenuItem doNeighborJoin = new JMenuItem ( "Infer Neighbor Joining Tree" );
                doNeighborJoin.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        if(e.getSource() == doNeighborJoin)
                        {

                            for(Proc1<String> listener : _neighborJoiningAnalysisRequestedListeners)
                            {
                                listener.execute(sequencingsTitleFinal);
                            }
                        }
                    }
                });
                analysisOptions.add (doNeighborJoin);

            }

        }
        else if(selections.length == 3)
        {
            final String firstEntryTitle =  (String) ((DefaultMutableTreeNode)selections[0].getLastPathComponent()).getUserObject();
            final boolean firstEntryTitleIsSequencingsData = _viewModel.isDataRecordTitleRepresentingSequencingsData(firstEntryTitle);
            boolean firstEntryTitleIsTraceData = _viewModel.isDataRecordTitleRepresentingTraceData(firstEntryTitle);

            final String secondEntryTitle = (String) ((DefaultMutableTreeNode)selections[1].getLastPathComponent()).getUserObject();
            boolean secondEntryTitleIsSequencingsData = _viewModel.isDataRecordTitleRepresentingSequencingsData(secondEntryTitle);
            boolean secondEntryTitleIsTraceData = _viewModel.isDataRecordTitleRepresentingTraceData(secondEntryTitle);

            final String thirdEntryTitle = (String) ((DefaultMutableTreeNode)selections[2].getLastPathComponent()).getUserObject();

            String sequencingsTitle = null, traceTitle = null, firstPositiveTitle = null;
            for(String title : Arrays.asList(firstEntryTitle, secondEntryTitle, thirdEntryTitle))
            {
                if(_viewModel.isDataRecordTitleRepresentingSequencingsData(title))
                {
                    sequencingsTitle = title;
                }
                else if(_viewModel.isDataRecordTitleRepresentingTraceData(title))
                {
                    traceTitle = title;
                }
                else if(_viewModel.isDataRecordTitleRepresentingFirstPositiveData(title))
                {
                    firstPositiveTitle = title;
                }
            }


            boolean showSnitkinAnalysisOption = sequencingsTitle != null && traceTitle != null && firstPositiveTitle != null;




           if(showSnitkinAnalysisOption)
           {
               final String sequencingsTitleFinal = sequencingsTitle, traceTitleFinal = traceTitle, firstPositiveTitleFinal = firstPositiveTitle;
               analysisOptions = new JPopupMenu ();
               final JMenuItem doSnitkin = new JMenuItem ( "Infer Snitkin Trans Map" );
               doSnitkin.addActionListener(new ActionListener() {
                   public void actionPerformed(ActionEvent e) {
                       if(e.getSource() == doSnitkin)
                       {

                           for(Proc3<String,String,String> listener : _snitkinTransMapAnalysisRequestedListeners)
                           {
                               listener.execute(sequencingsTitleFinal, traceTitleFinal, firstPositiveTitleFinal);
                           }
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

    private void onDataRecordAdded(WorkspaceVM.DataRecord record)
    {
        DefaultMutableTreeNode newEntryNode = new DefaultMutableTreeNode(record.Title);
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
                for(Proc1<File> addDataListener : _addDataListeners)
                {
                    addDataListener.execute(file);
                }
            }

        }
    }

}

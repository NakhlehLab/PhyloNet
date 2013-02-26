package edu.rice.cs.bioinfo.library.math.discrete;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/10/13
 * Time: 11:20 AM
 * To change this template use File | Settings | File Templates.
 */
public class Configurations<T> implements Iterable<List<T>>
{
    private class BinChoice
    {
        private final List<T> _options;

        private int _currentSelectionIndex = 0;

        private BinChoice(List<T> options) {
            _options = options;
        }

        public boolean canAdvance()
        {
            return _currentSelectionIndex < _options.size() - 1;
        }

        public void advance()
        {
            if(canAdvance())
            {
                _currentSelectionIndex++;
            }
            else
            {
                throw new IllegalStateException("Cannot advance bin choice because last option is selected.");
            }

        }

        public void reset()
        {
            _currentSelectionIndex = 0;
        }

        public T getChoice()
        {
            return _options.get(_currentSelectionIndex);
        }
    }

    private final List<List<T>> _bins;

    public Configurations(List<List<T>> bins)
    {
        _bins = bins;
    }

    public Iterator<List<T>> iterator()
    {
        final LinkedList<BinChoice> configuration = new LinkedList<BinChoice>();

        for(List<T> bin : _bins)
        {
            configuration.add(new BinChoice(bin));
        }

        return new Iterator<List<T>>()
        {
            private boolean  _hasNext = true;

            public boolean hasNext() {
                return _hasNext;
            }

            public List<T> next()
            {
                int numResets = 0;

                for (BinChoice bin : configuration)
                {
                    if(bin.canAdvance())
                    {
                        bin.advance();
                        return makeList(configuration);
                    }
                    else
                    {
                        bin.reset();
                        numResets++;
                    }
                }

                _hasNext = false;
                return makeList(configuration);

            }

            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    private LinkedList<T> makeList(LinkedList<BinChoice> configuration) {
        LinkedList<T> configList = new LinkedList<T>();

        for(BinChoice bin : configuration)
        {
            configList.add(bin.getChoice());
        }

        return configList;
    }
}

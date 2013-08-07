package containers;

public class IntArrayFactory implements IArrayFactory<Integer> {

    @Override
    public IArray<Integer> make(final int length) {
        return new IArray<Integer>() {

            private int[] array = new int[length];

            @Override
            public void set(int index, Integer type) {
                array[index] = type;
            }

            @Override
            public int getLength() {
                return length;
            }

            @Override
            public Integer get(int index) {
                return array[index];
            }
        };
    }
}

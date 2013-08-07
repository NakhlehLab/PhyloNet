//package containers;
//
//import java.util.ArrayList;
//
//public class ListOfArrayArrayFactory<Type> implements IArrayFactory<IListOfArrays<Type>> {
//	@Override
//	public IArray<IListOfArrays<Type>> make(final int length) {
//		return new IArray<IListOfArrays<Type>>() {
//			ArrayList<IListOfArrays<Type>> list = new ArrayList<IListOfArrays<Type>>(length);
//
//			@Override
//			public IListOfArrays<Type> get(int index) {
//				return list.get(index);
//			}
//
//			@Override
//			public void set(int index, IListOfArrays<Type> type) {
//				list.set(index, type);
//			}
//
//			@Override
//			public int getLength() {
//				return length;
//			}
//		};
//	}
//}

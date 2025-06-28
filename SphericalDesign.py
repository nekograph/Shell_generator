import itertools
import math
from collections import Counter


class E_n:
    
    def __init__(self, norm_sq: int):
        """
        E8格子の指定されたノルム二乗nのシェル(E_8)_nに属するベクトルをすべて返します。
        """
        if not isinstance(norm_sq, int) or norm_sq < 0:
            raise ValueError("エラー: norm_sqは非負整数である必要があります。")
        
        self.norm_sq = norm_sq 
        self._vectors = None 
        self._result_set = set() # E_8のシェルを入れる

    def get_shell(self) -> list[tuple]:
        """
        シェルを生成し、ソートされたベクトルのリストを返す。
        """
        if self._vectors is None: 
            self._generate_all_vectors() 
            self._vectors = sorted(list(self._result_set)) 
        return self._vectors 

    def _generate_all_vectors(self):
        """整数と半整数の両方のベクトル生成プロセスを開始する。"""
        # 1. 整数ベクトルの探索
        max_int_val = int(math.sqrt(self.norm_sq))
        self._find_integer_partitions(self.norm_sq, 8, max_int_val, [])

        # 2. 半整数ベクトルの探索
        if self.norm_sq > 0:
            target_4n = 4 * self.norm_sq
            max_odd_val = int(math.sqrt(target_4n))
            self._find_half_integer_partitions(target_4n, 8, max_odd_val, [])
    
    def _find_integer_partitions(self, target, k, max_val, partition):
        """ nを最大8個の平方数の和に分割する組み合わせを探す。"""
        if target == 0:
            # 組み合わせが決まったら_generate_vectors_from_partitionで並び替えや符号の組み合わせを探す
            self._generate_vectors_from_partition(partition + [0] * (8 - len(partition)), is_integer=True) # 余ったら０を入れる
            return
        if k == 0 or target < 0: # 見つからなかった
            return
        
        
        limit = int(math.sqrt(target))
        search_max = min(max_val, limit)
        for val in range(search_max, 0, -1):
            # 一つ成分が決まったら，ノルムをその成分分ひいて探す
            self._find_integer_partitions(target - val*val, k - 1, val, partition + [val])

    def _find_half_integer_partitions(self, target_4n, k, max_odd_val, partition):
        """ 4nを8個の奇数の平方数の和に分割する組み合わせを探す。"""
        if k == 0:
            if target_4n == 0:
                self._generate_vectors_from_partition(partition, is_integer=False)
            return
        if target_4n < 0:
            return

        limit = int(math.sqrt(target_4n))
        search_max = min(max_odd_val, limit)
        # 奇数のみを探すのでsearch_maxを超えない奇数から始めるようにする
        if search_max % 2 == 0:
            search_max -= 1
        for odd_val in range(search_max, 0, -2):
            self._find_half_integer_partitions(target_4n - odd_val*odd_val, k - 1, odd_val, partition + [odd_val])

    def _generate_vectors_from_partition(self, abs_coords, is_integer):
        """分割（絶対値リスト）から順列と符号の組み合わせを生成する。"""
        for p in set(itertools.permutations(abs_coords)):
            non_zero_indices = [i for i, v in enumerate(p) if v != 0] # 成分が0でない添え字を抽出

            # 各０でない成分に\pm1倍したものもvec_listに入れる
            for signs in itertools.product([-1, 1], repeat=len(non_zero_indices)):
                vec_list = list(p) # 計算可能にするためにリストに
                for i, sign in enumerate(signs):
                    vec_list[non_zero_indices[i]] *= sign 
                
                vec = tuple(vec_list)

                # 整数ベクトルと半整数ベクトルのE_8の元になる条件
                if is_integer:
                    if sum(vec) % 2 == 0: # 成分の和が２の倍数
                        self._result_set.add(vec)
                else: # 半整数
                    if sum(vec) % 4 == 0: # 2*xの成分の和が4の倍数
                        self._result_set.add(tuple(v / 2.0 for v in vec))
    

    # --- Pythonicな追加機能 ---
    def __len__(self):
        """len(instance)でベクトルの数を返せるようにする。"""
        return len(self.get_shell())

    def __iter__(self):
        """for vec in instance のようにイテレーションを可能にする。"""
        return iter(self.get_shell())

    def __getitem__(self, index):
        """instance[i] のようにインデックスでベクトルにアクセスできるようにする。"""
        return self.get_shell()[index]

# 整数ベクトルのサブクラス
class IntE_n(E_n):
    
    
    def get_shell(self) -> list[tuple]:
        if self._vectors is None:
            #  整数ベクトルの探索処理のみを呼び出す 
            max_int_val = int(math.sqrt(self.norm_sq))
            self._find_integer_partitions(self.norm_sq, 8, max_int_val, [])
            
            self._vectors = sorted(list(self._result_set))
        return self._vectors


# 半整数ベクトル用のサブクラス
class HalfE_n(E_n):
    

       # 親クラスのget_shellを上書きする
    def get_shell(self) -> list[tuple]:
        if self._vectors is None:
            #  半整数ベクトルの探索処理のみを呼び出す 
            if self.norm_sq > 0:
                target_4n = 4 * self.norm_sq
                max_odd_val = int(math.sqrt(target_4n))
                self._find_half_integer_partitions(target_4n, 8, max_odd_val, [])
            
            self._vectors = sorted(list(self._result_set))
        return self._vectors

# nの約数３乗和を求める関数
def sigma_3(n):
    if not isinstance(n,int) or n<0:
       raise ValueError("正整数を入力してください")
    
    sum = 0
    for i in range(int(n)):
        if n % (i+1) == 0:
            sum += (i+1)**3

    return sum

#-----------------
# 実行部分

if __name__ == '__main__':
    
    for n in range(2,20,2):
        gen = E_n(n)
        gen_int = IntE_n(n)
        gen_half = HalfE_n(n)
        formula = 16*(sigma_3(n)-sigma_3(int(n/2))) # 半整数の数を求める公式
        # num = 0
        # num_2 = 0
        # harm_positive = 0
        # harm_negative = 0
        # for i in range(len(gen_half)):
        #     gen_half_int = [int(2*gen_half[i][j]) for j in range(8)]
        #     harm_half = math.prod(gen_half_int[j] for j in range(8))
        #     if harm_half > 0:
        #         harm_positive += harm_half
        #         # print(f"正{gen_half_int},{harm_half}")
        #     else:
        #         harm_negative += harm_half
        #         # print(f"負{gen_half_int},{harm_half}")
        
        # print(len(gen_half),formula)
      

        print(n)
        # 整数ベクトルでharmを求める
        harm_1=0
        for i in range(len(gen_int)):
            harm_1 += math.prod(gen_int[i][j] for j in range(8)) #成分を全てかける関数
        print(f"整数:{harm_1}, {len(gen_int)}") 
        
        # 半整数ベクトルでharmを求める
        harm_2=0
        for i in range(len(gen_half)):
            harm_2 += math.prod(gen_half[i][j] for j in range(8))
        print(f"半整数:{harm_2},{len(gen_half)}")

        # 合計
        harm = 0
        for i in range(len(gen)):
            harm += math.prod(gen[i][j] for j in range(8))
        print(f"合計:{harm}, {len(gen)}")


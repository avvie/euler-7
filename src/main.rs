use std::{cmp::min, collections::HashMap, ops::{Index, IndexMut}, sync::{LazyLock, Mutex, atomic::{AtomicI64, Ordering}}, time::Instant};
use num_traits::pow::Pow;

static KEYS_FI_ARRAY: LazyLock<Mutex<HashMap<i64, Vec<i64>>>> = LazyLock::new(|| {
    Mutex::new(HashMap::new())
});

pub struct FIArray{
    x : i64,
    isqrt: i64,
    arr: Vec<i64>
}

impl FIArray {
    pub fn new(x:i64) -> Self{
        let sqrt = isqrt(x);
        let mut l = 2* sqrt;
        if sqrt == floor_division(x, sqrt) {
            l -=1;
        }

        FIArray { x, isqrt: sqrt, arr: vec![0;l as usize] }
    }

    pub fn index_of(self, v :i64) -> i64 {
        if v <= self.isqrt {
            return v-1;
        }
        return self.arr.len() as i64 - floor_division(self.x, v);
    }
}

impl IndexMut<i64> for FIArray {
    fn index_mut(&mut self, v: i64) -> &mut Self::Output {
        if v <= self.isqrt {
            &mut self.arr[(v - 1) as usize]
        } else {
            let idx = (self.x / v) as usize;
            let len = self.arr.len();
            &mut self.arr[len- idx]
        }
    }
}

impl Index<i64> for FIArray {
    type Output = i64;

    fn index(&self, v: i64) -> &Self::Output {
        if v <= 0 {
            return &0;
        }
        if v <= self.isqrt {
            &self.arr[(v - 1) as usize]
        } else {
            let idx = (self.x / v) as usize;
            &self.arr[self.arr.len() - idx]
        }
    }
}

fn keys_fi(x: i64) -> Vec<i64> {
    if let Some(v) = KEYS_FI_ARRAY.lock().unwrap().get(&x) {
        return v.clone();
    }

    let rt = isqrt(x);
    let mut result = Vec::new();

    // Add 1 to rt
    for j in 1..=rt {
        result.push(j);
    }

    // Add x // j for j from rt down to 1, but skip duplicate at boundary
    if rt != x / rt {
        result.push(x / rt);
    }

    for j in (1..rt).rev() {
        result.push(x / j);
    }

    KEYS_FI_ARRAY.lock().unwrap().insert(x, result.clone());
    result
}   

#[derive(Debug, Clone)]
pub struct Fenwick<T> {
    arr: Vec<T>,
}

impl<T> Fenwick<T>
where
    T: Copy + std::ops::Add<Output = T> + std::ops::Sub<Output = T> + Default
    + std::ops::AddAssign,
{
    pub fn new(len: usize, def : T) -> Self {
        let mut arr = vec![def; len];
        for i in 1..len{
            let j = i + (i & (!i +1));
            if j <= len{
                let val = arr[i-1];
                arr[j-1] += val;
            }
        }
        Self {
            arr: arr,
        }
    }

    /// Returns the length of the underlying array.
    pub fn len(&self) -> usize {
        self.arr.len()
    }

    /// Computes the prefix sum from index 0 to i (inclusive).
    /// Time complexity: O(log i)
    pub fn sum(&self, mut i: usize) -> T {
        let mut sum = T::default();
        i += 1; // Convert to 1-based indexing
        while i > 0 {
            sum = sum + self.arr[i - 1]; // Access 0-based array
            i -= i & (!i + 1); // Clear the lowest set bit
        }
        sum
    }

    /// Adds `delta` to the element at index `i`.
    /// Time complexity: O(log i)
    pub fn add_to(&mut self, mut i: usize, delta: T) {
        i += 1; // Convert to 1-based indexing
        while i <= self.arr.len() {
            self.arr[i - 1] = self.arr[i - 1] + delta;
            i += i & (!i + 1); // Move to next index with lowest set bit added
        }
    }

    /// Gets the value at index `i` of the original array.
    /// Time complexity: O(log i)
    pub fn get(&self, i: usize) -> T {
        if i == 0 {
            self.sum(0)
        } else {
            self.sum(i) - self.sum(i - 1)
        }
    }

    /// Sets the value at index `i` to `value`.
    /// Time complexity: O(log i)
    pub fn set(&mut self, i: usize, value: T) {
        let current = self.get(i);
        self.add_to(i, value - current);
    }
}

fn isqrt(x:i64) -> i64 {
    return (f64::sqrt(x as f64)).floor() as i64;
}

fn floor_division(x:i64, n: i64)-> i64 {
    return ((x as f64) / (n as f64)).floor() as i64;
}


pub fn lucy_fenwich(x:i64, f: fn(i64)) -> FIArray {
    let mut s = FIArray::new(x);

    let xf: f64 = x as f64;
    let mut y :i64;

    if x==1 {
        y = 1;
    }
    else{
        y = num_traits::float::FloatCore::round(0.35 * (xf.pow(2.0/3.0) / xf.ln().pow(2.0 /3.0)) ) as i64;
        y = y.min(4e9 as i64);
        y = y.max(s.isqrt+1);
    }
    let mut sieve_raw = vec![false; (y+1) as usize];
    let mut sieve = Fenwick::<i64>::new((y+1) as usize, 0);
    sieve_raw[0] = true;
    sieve_raw[1] = true;
    sieve.set(1, 0);
    sieve.set(0, 0);
    let x_div_y = floor_division(x, y);

    for i in 2..=y{
        sieve.add_to(i as usize, 1);
    }


    let v_h = keys_fi(x);
    for (i, &v) in v_h.iter().enumerate() {
        if v > 1 {
            s.arr[i] = v-1;
        }
        else{
            s.arr[i] = 0;
        }
    }

    for p in 2..=s.isqrt {
        if !sieve_raw[p as usize] {
            f(p);
            let sp = sieve.sum((p-1) as usize);
            let j = p*p;
            let lim = min(x_div_y, floor_division(x, j));

            let len = s.arr.len();
            s.arr[len-1] -= s0(floor_division(x, p), y, &sieve, &s) -sp;
            for i in p..lim {
                if sieve_raw[i as usize] {
                    continue;
                }
                let count = s0(floor_division(x, i*p), y,  &sieve, &s) -sp;
                f(p*count);
                s.arr[len-i as usize] -= count;
            }

            let mut j_index = j;
            while  j_index <= y {
                if !sieve_raw[j_index as usize]{
                    sieve_raw[j_index as usize] = true;
                    f(p);
                    sieve.add_to(j_index as usize, -1);
                }
                j_index += p;
            }
        }
    }
    for (i, &v) in v_h.iter().enumerate() {
        if v > y {
            break;
        }
        s.arr[i] = sieve.sum(v as usize);
    }
    println!("sieve sum y: {}", sieve.sum(y as usize));
    return s;
}

fn s0(x:i64, y: i64, fen: &Fenwick<i64>, s: &FIArray) -> i64 {
    if x <= y {
        return fen.sum(x as usize);
    }
    return s[x];
}

static SUM: AtomicI64 = AtomicI64::new(0);
fn add_mod_9(x: i64){
    SUM.fetch_add(x, Ordering::Relaxed);

}

fn main() {
    println!("floor div check: {}", floor_division(10,2));
    let n = 1000000000000;
    let start = Instant::now();
    let s = lucy_fenwich(n, add_mod_9);
    let elapsed = start.elapsed();
    println!("Elapsed time: {:?}", elapsed);
    println!("the sum mod9 is {}", {SUM.load(Ordering::Relaxed)});
    println!("Lucy counted {}", s[n]);

    // println!("THere are {} primes in range.",primal::StreamingSieve::prime_pi(n as usize));
    // let elapsed = start.elapsed();
    // println!("Elapsed time: {:?}", elapsed);

}
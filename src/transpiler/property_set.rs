//! `PropertySet`: keyed, type-erased metadata passed between passes.

use std::any::Any;
use std::collections::HashMap;
use std::fmt;

/// A type-safe heterogeneous map for transpiler-pass metadata.
///
/// Keys are strings; values are downcast on retrieval. Passes typically agree
/// on property names as an out-of-band contract.
#[derive(Default)]
pub struct PropertySet {
    properties: HashMap<String, Box<dyn Any>>,
}

impl fmt::Debug for PropertySet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let keys: Vec<&str> = self.properties.keys().map(String::as_str).collect();
        f.debug_struct("PropertySet")
            .field("keys", &keys)
            .finish()
    }
}

impl PropertySet {
    /// Creates an empty [`PropertySet`].
    pub fn new() -> Self {
        Self::default()
    }

    /// Inserts a typed value, replacing any previous value at `key`.
    pub fn insert<T: 'static>(&mut self, key: &str, value: T) {
        self.properties.insert(key.to_string(), Box::new(value));
    }

    /// Retrieves a reference to a value of type `T`, returning `None` if the
    /// key is missing or the stored type does not match.
    pub fn get<T: 'static>(&self, key: &str) -> Option<&T> {
        self.properties
            .get(key)
            .and_then(|any| any.downcast_ref::<T>())
    }

    /// Removes and returns a value of type `T`. If the type does not match,
    /// the value is re-inserted and `None` is returned.
    pub fn remove<T: 'static>(&mut self, key: &str) -> Option<Box<T>> {
        if let Some(any) = self.properties.remove(key) {
            match any.downcast::<T>() {
                Ok(v) => return Some(v),
                Err(orig) => {
                    self.properties.insert(key.to_string(), orig);
                }
            }
        }
        None
    }

    /// Clears all properties.
    pub fn clear(&mut self) {
        self.properties.clear();
    }

    /// True iff the set contains `key`.
    pub fn contains(&self, key: &str) -> bool {
        self.properties.contains_key(key)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_property_set_insertion_and_retrieval() {
        let mut p = PropertySet::new();
        p.insert("depth", 15_usize);
        p.insert("name", "c1".to_string());
        assert_eq!(p.get::<usize>("depth"), Some(&15));
        assert_eq!(p.get::<String>("name"), Some(&"c1".to_string()));
        assert_eq!(p.get::<bool>("missing"), None);
        assert_eq!(p.get::<usize>("name"), None);
    }

    #[test]
    fn test_property_set_remove_type_mismatch() {
        let mut p = PropertySet::new();
        p.insert("x", 42_usize);
        assert!(p.remove::<String>("x").is_none());
        assert_eq!(p.get::<usize>("x"), Some(&42));
    }
}
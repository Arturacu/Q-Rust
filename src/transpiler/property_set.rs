use std::any::Any;
use std::collections::HashMap;

/// A generalized metadata storage for passing information between transpiler passes.
/// This prevents expensive recalculations (e.g., rescanning circuits for properties).
#[derive(Default)]
pub struct PropertySet {
    properties: HashMap<String, Box<dyn Any>>,
}

impl PropertySet {
    /// Creates a new, empty PropertySet.
    pub fn new() -> Self {
        Self {
            properties: HashMap::new(),
        }
    }

    /// Inserts a property into the set.
    pub fn insert<T: 'static>(&mut self, key: &str, value: T) {
        self.properties.insert(key.to_string(), Box::new(value));
    }

    /// Retrieves an immutable reference to a property if it exists and matches the typed signature.
    pub fn get<T: 'static>(&self, key: &str) -> Option<&T> {
        self.properties
            .get(key)
            .and_then(|any| any.downcast_ref::<T>())
    }

    /// Removes and retrieves a property from the set.
    pub fn remove<T: 'static>(&mut self, key: &str) -> Option<Box<T>> {
        if let Some(any) = self.properties.remove(key) {
            if any.is::<T>() {
                return any.downcast::<T>().ok();
            }
            // If type doesn't match, carefully re-insert it rather than dropping silently
            self.properties.insert(key.to_string(), any);
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_property_set_insertion_and_retrieval() {
        let mut props = PropertySet::new();
        props.insert("depth", 15_usize);
        props.insert("name", "circuit_1".to_string());

        assert_eq!(props.get::<usize>("depth"), Some(&15));
        assert_eq!(props.get::<String>("name"), Some(&"circuit_1".to_string()));
        assert_eq!(props.get::<bool>("missing"), None);
        assert_eq!(props.get::<usize>("name"), None); // Type mismatch
    }
}

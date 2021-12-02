use regex::Regex;
use std::error::Error;
use std::fmt;

#[derive(Debug, PartialEq, Clone)]
pub struct OptionalRegion {
    pub path: String,
    pub start: Option<u64>,
    pub end: Option<u64>,
}

impl fmt::Display for OptionalRegion {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Use `self.number` to refer to each positional data point.
        match self.start {
            Some(start) => match self.end {
                Some(end) => write!(f, "{}:{}-{}", self.path, start, end),
                None => write!(f, "{}:{}", self.path, start),
            },
            None => write!(f, "{}", self.path),
        }
    }
}

impl OptionalRegion {
    pub fn interval(&self) -> Option<u64> {
        if let Some(start) = self.start {
            if let Some(end) = self.end {
                if start < end {
                    return Some(end - start);
                } else {
                    return Some(start - end);
                }
            }
        }
        None
    }

    pub fn inverted(&self) -> Option<bool> {
        if let Some(_start) = self.start {
            if let Some(_end) = self.end {
                return Some(self.start > self.end);
            }
        }
        None
    }

    pub fn new_with_prefix(path: String, chr_prefix: &str) -> Result<Self, Box<dyn Error>> {
        let re = Regex::new(r"^(.+):(\d*)-?(\d*)$").unwrap();
        let caps = re.captures(&path).ok_or("Invalid genomic range")?;
        let mut path_str = caps.get(1).ok_or("Parse Path Error")?.as_str();

        let path_string: String;
        if chr_prefix.len() == 0 {
            if path_str.starts_with("chr") {
                path_str = &path_str[3..];
            }
            path_string = path_str.to_string();
        } else {
            if path_str.starts_with(chr_prefix) {
                path_str = &path_str[chr_prefix.len()..];
            }
            if path_str.len() < chr_prefix.len() {
                path_string = format!("{}{}", chr_prefix, path_str);
            } else {
                path_string = path_str.to_string()
            }
        }
        let start = caps.get(2).and_then(|t| t.as_str().parse::<u64>().ok());
        let end = caps.get(3).and_then(|t| t.as_str().parse::<u64>().ok());
        return Ok(OptionalRegion {
            path: path_string,
            start: start,
            end: end,
        });
    }

    pub fn new(path: &str) -> Result<Self, Box<dyn Error>> {
        let re = Regex::new(r"^(.+):(\d*)-?(\d*)$").unwrap();
        let caps = re.captures(path).ok_or("Invalid genomic range")?;
        let path = caps.get(1).ok_or("Parse Path Error")?;
        let start = caps.get(2).and_then(|t| t.as_str().parse::<u64>().ok());
        let end = caps.get(3).and_then(|t| t.as_str().parse::<u64>().ok());
        return Ok(OptionalRegion {
            path: path.as_str().to_string(),
            start: start,
            end: end,
        });
    }

    pub fn uuid(self: &OptionalRegion) -> String {
        return format!("{}", self);
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct StringRegion {
    pub path: String, // Requires no prefix
    pub start: u64,
    pub end: u64,
    inverted: bool,
}

impl fmt::Display for StringRegion {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Use `self.number` to refer to each positional data point.
        if self.inverted {
            write!(f, "{}:{}-{}", self.path, self.end, self.start)
        } else {
            write!(f, "{}:{}-{}", self.path, self.start, self.end)
        }
    }
}

impl StringRegion {
    pub fn interval(&self) -> u64 {
        return self.end - self.start;
    }
    pub fn inverted(&self) -> bool {
        self.inverted
    }
    pub fn start(&self) -> u64 {
        self.start
    }
    pub fn left(&self) -> u64 {
        if self.inverted {
            return self.end
        } else {
            return self.start
        }
    }
    pub fn right(&self) -> u64 {
        if self.inverted {
            return self.start
        } else {
            return self.end
        }
    }

    pub fn end(&self) -> u64 {
        self.end
    }

    pub fn extend(&mut self, len: u64) {
        self.start -= len;
        if self.start <= 0 {
            self.start = 0;
        }
        self.end += len;
    }
    // It is used on converting dna-sequence region to bed-style region.
    pub fn start_minus(&mut self) {
        self.start = self.start - 1;
    }

    pub fn new_with_prefix(path: String, chr_prefix: &str) -> Result<Self, Box<dyn Error>> {
        let re = Regex::new(r"^(.+):(\d+)-?(\d*)$").unwrap();
        let caps = re.captures(&path).ok_or("Invalid genomic range")?;
        let mut path_str = caps.get(1).ok_or("Parse Path Error")?.as_str();
        let path_string: String;
        if chr_prefix.len() == 0 {
            if path_str.starts_with("chr") {
                path_str = &path_str[3..];
            }
            path_string = path_str.to_string();
        } else {
            if path_str.starts_with(chr_prefix) {
                path_str = &path_str[chr_prefix.len()..]; // .replace("chr", "");
            }
            if path_str.len() < chr_prefix.len() {
                path_string = format!("{}{}", chr_prefix, path_str);
            } else {
                path_string = path_str.to_string()
            }
        }
        let start = caps.get(2).ok_or("Parse Start Position Error")?;
        let end = caps.get(3).ok_or("Parse end Position Error")?;
        let start_str: &str = start.as_str().as_ref();
        let end_str: &str = end.as_str().as_ref();
        let start_u64: u64 = start_str
            .parse::<u64>()
            .map_err(|e| "Parse Int Error, ".to_string() + &e.to_string())?;
        let end_u64: u64 = end_str
            .parse::<u64>()
            .map_err(|e| "Parse Int Error, ".to_string() + &e.to_string())?;
        Ok(StringRegion::new_inner(path.to_string(), start_u64, end_u64))
    }

    fn new_regexp(path: &str) -> Result<Self, Box<dyn Error>> {
        let re = Regex::new(r"^(.+):(\d+)-?(\d*)$").unwrap();
        let caps = re.captures(path).ok_or("Invalid genomic range")?;
        let path = caps.get(1).ok_or("Parse Path Error")?;
        let start = caps.get(2).ok_or("Parse Start Position Error")?;
        let end = caps.get(3).ok_or("Parse end Position Error")?;
        let start_str: &str = start.as_str().as_ref();
        let end_str: &str = end.as_str().as_ref();
        let start_u64: u64 = start_str
            .parse::<u64>()
            .map_err(|e| "Parse Int Error, ".to_string() + &e.to_string())?;
        let end_u64: u64 = end_str
            .parse::<u64>()
            .map_err(|e| "Parse Int Error, ".to_string() + &e.to_string())?;
        Ok(StringRegion::new_inner(path.as_str().to_string(), start_u64, end_u64))
    }

    pub fn new(path: &str) -> Result<Self, Box<dyn Error>> {
        let caps: Vec<&str> = path.split_whitespace().collect();
        if caps.len() < 3 {
            return StringRegion::new_regexp(path);
        }
        let path = caps[0];
        let start = caps[1];
        let end = caps[2];
        let start_u64: u64 = start
            .parse::<u64>()
            .map_err(|e| "Parse Int Error, ".to_string() + &e.to_string())?;
        let end_u64: u64 = end
            .parse::<u64>()
            .map_err(|e| "Parse Int Error, ".to_string() + &e.to_string())?;
        Ok(StringRegion::new_inner(path.to_string(), start_u64, end_u64))
    }

    pub fn new_inner(path: String, start_u64: u64, end_u64: u64) -> Self {
        if start_u64 > end_u64 {
            StringRegion {
                path: path.to_string(),
                start: end_u64,
                end: start_u64,
                inverted: true
            }
        } else {
            StringRegion {
                path: path.to_string(),
                start: start_u64,
                end: end_u64,
                inverted: false
            }
        }
    }

    pub fn uuid(&self) -> String {
        return format!("{}", self);
    }
}

#[derive(Clone, PartialEq, Debug)]
pub struct Region {
    ref_id: u64,
    start: u64,
    end: u64,
}

impl Region {
    /// Creates new region. `ref_id` is 0-based, `start-end` is 0-based half-open interval.
    pub fn new(ref_id: u64, start: u64, end: u64) -> Region {
        assert!(
            start <= end,
            "Region: start should not be greater than end ({} > {})",
            start,
            end
        );
        Region { ref_id, start, end }
    }

    pub fn convert<F>(
        path: &StringRegion,
        to_id: F,
    ) -> std::result::Result<Self, Box<dyn std::error::Error>>
    where
        F: Fn(&str) -> Option<u64>,
    {
        Ok(Region {
            ref_id: to_id(&path.path).ok_or("Error: the reference id is not recognized.")?,
            start: path.start,
            end: path.end,
        })
    }

    pub fn parse<F>(path: &str, to_id: F) -> std::result::Result<Self, Box<dyn std::error::Error>>
    where
        F: Fn(&str) -> Option<u64>,
    {
        let re = Regex::new(r"^(.+):(\d*)-?(\d*)$").unwrap();
        let caps = re.captures(path).ok_or("Invalid genomic range")?;
        let path = caps
            .get(1)
            .and_then(|t| Some(t.as_str()))
            .ok_or("Parse Path Error")?;
        let start = caps
            .get(2)
            .and_then(|t| t.as_str().parse::<u64>().ok())
            .ok_or("Error: the reference start is not recognized.")?;
        let end = caps
            .get(3)
            .and_then(|t| t.as_str().parse::<u64>().ok())
            .ok_or("Error: the reference end is not recognized.")?;

        return Ok(Region {
            ref_id: to_id(path).ok_or("Error: the reference id is not recognized.")?,
            start: start,
            end: end,
        });
    }

    pub fn ref_id(&self) -> u64 {
        self.ref_id
    }

    pub fn start(&self) -> u64 {
        self.start
    }

    pub fn end(&self) -> u64 {
        self.end
    }

    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    pub fn set_ref_id(&mut self, ref_id: u64) {
        self.ref_id = ref_id;
    }

    pub fn set_start(&mut self, start: u64) {
        assert!(
            start <= self.end,
            "Region: start should not be greater than end ({} > {})",
            start,
            self.end
        );
        self.start = start;
    }

    pub fn set_end(&mut self, end: u64) {
        assert!(
            self.start <= end,
            "Region: start should not be greater than end ({} > {})",
            self.start,
            end
        );
        self.end = end;
    }

    pub fn contains(&self, ref_id: u64, pos: u64) -> bool {
        self.ref_id == ref_id && self.start <= pos && pos < self.end
    }

    pub fn include(&self, range: &Region) -> bool {
        self.ref_id == range.ref_id && self.start <= range.start && range.end < self.end
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region_format(path: &str) -> String {
        format!("{}", StringRegion::new(path).unwrap())
    }

    #[test]
    fn region_works() {
        assert_eq!(StringRegion::new("").ok(), None);
        assert_eq!(StringRegion::new(":10-20").ok(), None);
        assert_eq!(
            StringRegion::new("chr1:12000-12001").ok(),
            Some(StringRegion {
                path: "chr1".to_string(),
                start: 12000,
                end: 12001
            })
        );
        assert_eq!(
            StringRegion::new("chr1:1200943-1201000").ok(),
            Some(StringRegion {
                path: "chr1".to_string(),
                start: 1200943,
                end: 1201000
            })
        );
    }

    #[test]
    fn region_format_works() {
        let a = "chr1:12000-12001";
        assert_eq!(region_format(a), a);
        let b = "10:120-120001";
        assert_eq!(region_format(b), b);
    }
}

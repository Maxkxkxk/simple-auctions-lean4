import SimpleAuctions.Kernel

/-! The two kernel evaluations (heavy: big-integer arithmetic inside the Lean kernel). -/
namespace SimpleAuctions
namespace Cert

set_option maxRecDepth 100000 in
theorem revCheck_true : revCheck = true := by decide +kernel

set_option maxRecDepth 100000 in
theorem brevCheck_true : brevCheck = true := by decide +kernel

end Cert
end SimpleAuctions

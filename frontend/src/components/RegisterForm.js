import React, { useState, useContext } from "react";
import { AuthContext } from "../AuthContext";
import { Icons } from "./Icons";

const RegisterForm = ({ toggleForm }) => {
  const [email, setEmail] = useState("");
  const [password, setPassword] = useState("");
  const [confirmPassword, setConfirmPassword] = useState("");
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState("");

  const { register } = useContext(AuthContext);

  const handleSubmit = async (e) => {
    e.preventDefault();
    setError("");
    
    if (!email || !password || !confirmPassword) {
      setError("All fields are required");
      return;
    }

    if (password !== confirmPassword) {
      setError("Passwords do not match");
      return;
    }

    if (password.length < 8) {
      setError("Password must be at least 8 characters long");
      return;
    }

    setIsLoading(true);

    try {
      const result = await register(email, password);
      if (!result.success) {
        setError(result.error);
      }
    } catch (err) {
      setError("An unexpected error occurred");
      console.error("Registration error:", err);
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <div className="register-form" style={{ padding: "30px" }}>
      <h3 style={{ 
        fontSize: "24px", 
        fontWeight: "700", 
        marginBottom: "20px", 
        display: "flex", 
        alignItems: "center", 
        gap: "10px" 
      }}>
        <Icons.verified /> Create Account
      </h3>

      {error && (
        <div style={{ 
          color: "#f43f5e", 
          marginBottom: "15px", 
          padding: "10px", 
          borderRadius: "6px",
          backgroundColor: "rgba(244, 63, 94, 0.1)",
          border: "1px solid rgba(244, 63, 94, 0.2)"
        }}>
          {error}
        </div>
      )}

      <form onSubmit={handleSubmit}>
        <div style={{ marginBottom: "15px" }}>
          <label 
            htmlFor="email" 
            style={{ 
              display: "block", 
              marginBottom: "5px", 
              fontSize: "14px", 
              fontWeight: "500", 
              color: "#94a3b8" 
            }}
          >
            Email Address
          </label>
          <div style={{ position: "relative" }}>
            <span style={{ 
              position: "absolute", 
              left: "12px", 
              top: "50%", 
              transform: "translateY(-50%)", 
              color: "#94a3b8" 
            }}>
              <Icons.info />
            </span>
            <input
              id="email"
              type="email"
              value={email}
              onChange={(e) => setEmail(e.target.value)}
              style={{ 
                width: "100%", 
                padding: "12px 12px 12px 36px", 
                backgroundColor: "rgba(15, 23, 42, 0.7)", 
                border: "1px solid rgba(255, 255, 255, 0.1)", 
                borderRadius: "6px", 
                color: "#f0f4f8", 
                fontSize: "14px" 
              }}
              placeholder="your@email.com"
              required
            />
          </div>
        </div>

        <div style={{ marginBottom: "15px" }}>
          <label 
            htmlFor="password" 
            style={{ 
              display: "block", 
              marginBottom: "5px", 
              fontSize: "14px", 
              fontWeight: "500", 
              color: "#94a3b8" 
            }}
          >
            Password
          </label>
          <div style={{ position: "relative" }}>
            <span style={{ 
              position: "absolute", 
              left: "12px", 
              top: "50%", 
              transform: "translateY(-50%)", 
              color: "#94a3b8" 
            }}>
              <Icons.lock />
            </span>
            <input
              id="password"
              type="password"
              value={password}
              onChange={(e) => setPassword(e.target.value)}
              style={{ 
                width: "100%", 
                padding: "12px 12px 12px 36px", 
                backgroundColor: "rgba(15, 23, 42, 0.7)", 
                border: "1px solid rgba(255, 255, 255, 0.1)", 
                borderRadius: "6px", 
                color: "#f0f4f8", 
                fontSize: "14px" 
              }}
              placeholder="••••••••"
              required
            />
          </div>
        </div>

        <div style={{ marginBottom: "20px" }}>
          <label 
            htmlFor="confirmPassword" 
            style={{ 
              display: "block", 
              marginBottom: "5px", 
              fontSize: "14px", 
              fontWeight: "500", 
              color: "#94a3b8" 
            }}
          >
            Confirm Password
          </label>
          <div style={{ position: "relative" }}>
            <span style={{ 
              position: "absolute", 
              left: "12px", 
              top: "50%", 
              transform: "translateY(-50%)", 
              color: "#94a3b8" 
            }}>
              <Icons.lock />
            </span>
            <input
              id="confirmPassword"
              type="password"
              value={confirmPassword}
              onChange={(e) => setConfirmPassword(e.target.value)}
              style={{ 
                width: "100%", 
                padding: "12px 12px 12px 36px", 
                backgroundColor: "rgba(15, 23, 42, 0.7)", 
                border: "1px solid rgba(255, 255, 255, 0.1)", 
                borderRadius: "6px", 
                color: "#f0f4f8", 
                fontSize: "14px" 
              }}
              placeholder="••••••••"
              required
            />
          </div>
        </div>

        <button
          type="submit"
          disabled={isLoading}
          style={{
            width: "100%",
            padding: "12px",
            borderRadius: "8px",
            backgroundColor: "rgba(16, 185, 129, 0.9)",
            color: "white",
            fontSize: "16px",
            fontWeight: "600",
            border: "none",
            cursor: isLoading ? "not-allowed" : "pointer",
            opacity: isLoading ? "0.7" : "1",
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
            gap: "8px",
            marginBottom: "20px"
          }}
        >
          {isLoading ? (
            <>
              <span className="spinner"></span>
              Creating Account...
            </>
          ) : (
            "Create Account"
          )}
        </button>

        <div style={{ textAlign: "center", fontSize: "14px", color: "#94a3b8" }}>
          Already have an account?{" "}
          <button
            type="button"
            onClick={toggleForm}
            style={{
              background: "none",
              border: "none",
              color: "#38bdf8",
              cursor: "pointer",
              fontWeight: "500",
              padding: "0",
              textDecoration: "underline"
            }}
          >
            Sign in
          </button>
        </div>
      </form>
    </div>
  );
};

export default RegisterForm;
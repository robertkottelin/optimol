# Environment Variables Configuration

This document describes all environment variables used in the OptimizeMolecule application.

## Backend Environment Variables

**Location**: `/root/optimol/backend/.env`

### Required Variables

#### JWT Configuration
- **`JWT_SECRET_KEY`** - Secret key for JWT token signing
  - **Type**: String (recommended: 64+ random characters)
  - **Example**: `zyb8Pgc8sg5pemJAt-CEaOPd_jiHgKnbtiRiz8y3Fp8xdj7i4v64eQEAW0kW-5YOhqYO2x2xSuVVmKiuFjKjQA`
  - **Where to get**: Generate using: `python -c "import secrets; print(secrets.token_urlsafe(64))"`

#### Stripe Configuration
- **`STRIPE_SECRET_KEY`** - Stripe API secret key (server-side)
  - **Type**: String (starts with `sk_live_` or `sk_test_`)
  - **Example**: `sk_live_51GhEaMJQZaUHxA2L...`
  - **Where to get**: [Stripe Dashboard → Developers → API Keys](https://dashboard.stripe.com/apikeys)
  - **⚠️ IMPORTANT**: Never expose this key to the client

- **`STRIPE_WEBHOOK_SECRET`** - Stripe webhook signing secret
  - **Type**: String (starts with `whsec_`)
  - **Example**: `whsec_GPDNifAEqIhUJMFp2M7Bkcgi18EGrhYz`
  - **Where to get**: [Stripe Dashboard → Developers → Webhooks](https://dashboard.stripe.com/webhooks)
  - **Note**: Create webhook endpoint at `https://yourdomain.com/api/webhook`

- **`STRIPE_SUBSCRIPTION_PRICE_ID`** - Stripe price ID for subscription
  - **Type**: String (starts with `price_`)
  - **Example**: `price_1R8QFRJQZaUHxA2Lgbz0WRd5`
  - **Where to get**: [Stripe Dashboard → Products](https://dashboard.stripe.com/products)
  - **Note**: Create a recurring price for your subscription product

### Testing with Stripe Test Mode

To use Stripe in test mode, replace the live keys with test keys:

```bash
# Test mode configuration
STRIPE_SECRET_KEY=sk_test_YOUR_TEST_KEY_HERE
STRIPE_WEBHOOK_SECRET=whsec_YOUR_TEST_WEBHOOK_SECRET_HERE
STRIPE_SUBSCRIPTION_PRICE_ID=price_YOUR_TEST_PRICE_ID_HERE
```

Test credit cards: [Stripe Testing Documentation](https://stripe.com/docs/testing)

---

## Frontend Environment Variables

**Location**: `/root/optimol/frontend/.env`

### Required Variables

- **`REACT_APP_STRIPE_PUBLISHABLE_KEY`** - Stripe publishable key (client-side)
  - **Type**: String (starts with `pk_live_` or `pk_test_`)
  - **Example**: `pk_live_Tfh90MeFSg6jVjRCaMExaGug0078PAfanh`
  - **Where to get**: [Stripe Dashboard → Developers → API Keys](https://dashboard.stripe.com/apikeys)
  - **Note**: Safe to expose publicly (client-side key)

### Optional Variables

- **`REACT_APP_API_BASE_URL`** - Backend API base URL
  - **Type**: String (URL path or full URL)
  - **Default**: `/api` (uses nginx proxy in production)
  - **Example for local dev**: `http://localhost:5000`
  - **Note**: Only needed for local development or if backend is on different domain

---

## Environment Variable Loading

### Backend (Flask)
- Loaded via systemd `EnvironmentFile` directive
- Configuration: `/etc/systemd/system/optimol-backend.service`
- After updating `.env`, restart: `sudo systemctl restart optimol-backend`

### Frontend (React)
- Loaded at build time by Create React App
- Variables must start with `REACT_APP_`
- After updating `.env`, rebuild: `npm run build`
- Deploy to nginx: `cp -r build/* /var/www/optimizemolecule/`

---

## Security Best Practices

1. **Never commit `.env` files to git** - They are gitignored by default
2. **Use `.env.example` files** - Provided as templates (no actual secrets)
3. **Rotate secrets regularly** - Especially JWT_SECRET_KEY
4. **Use test keys for development** - Switch to live keys only for production
5. **Limit webhook IP access** - Configure firewall to only accept webhooks from Stripe IPs

---

## Example Files

### Backend `.env.example`
See: `/root/optimol/backend/.env.example`

### Frontend `.env.example`
See: `/root/optimol/frontend/.env.example`

---

## Troubleshooting

### Backend not connecting to Stripe
- Check: `journalctl -u optimol-backend -n 50`
- Verify: `STRIPE_SECRET_KEY` is set and starts with `sk_`
- Test: `curl http://localhost:5000/health`

### Frontend Stripe not loading
- Check browser console for errors
- Verify: `REACT_APP_STRIPE_PUBLISHABLE_KEY` in `.env`
- Ensure: Frontend was rebuilt after changing `.env`
- Test: View page source and search for "pk_live" in JavaScript bundle

### Webhooks not working
- Check: `STRIPE_WEBHOOK_SECRET` is configured
- Verify: Webhook endpoint is `https://yourdomain.com/api/webhook`
- Test: Send test webhook from Stripe Dashboard

---

## Quick Reference: Where Keys Are Used

| Variable | File | Line(s) | Usage |
|----------|------|---------|-------|
| `STRIPE_SECRET_KEY` | `backend/user.py` | 18 | Initialize Stripe API |
| `STRIPE_SUBSCRIPTION_PRICE_ID` | `backend/user.py` | 21, 264, 390 | Create subscriptions |
| `STRIPE_WEBHOOK_SECRET` | `backend/user.py` | 536 | Verify webhook signatures |
| `JWT_SECRET_KEY` | `backend/app.py` | 36 | JWT token signing |
| `REACT_APP_STRIPE_PUBLISHABLE_KEY` | `frontend/src/App.js` | 32 | Initialize Stripe.js |
| `REACT_APP_API_BASE_URL` | `frontend/src/AuthContext.js` | 50 | API endpoint base URL |

---

Last Updated: October 27, 2025
